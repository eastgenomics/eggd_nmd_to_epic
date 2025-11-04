import dxpy
import concurrent.futures
import json
from collections import Counter

# Define current project as current workspace
current_project_id = dxpy.WORKSPACE_ID
current_project = dxpy.api.project_describe(current_project_id)
projects = [current_project]

# Class for NMD processing
class NMDProcessor:
    @staticmethod
    def call_in_parallel(func, items, ignore_missing=False, **kwargs):
        """
        Calls the given function in parallel using concurrent.futures on
        the given set of items (i.e for calling dxpy.describe() on multiple
        object IDs).

        Additional arguments specified to kwargs are directly passed to the
        specified function.

        Parameters
        ----------
        func : callable
            function to call on each item
        items : list
            list of items to call function on
        ignore_missing : bool
            controls if to just print a warning instead of raising an
            exception on a dxpy.exceptions.ResourceNotFound being raised.
            This is most likely from a file that has been deleted and we are
            just going to default to ignoring these

        Returns
        -------
        list
            list of responses
        """
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
            futures = {executor.submit(func, item, **kwargs): item for item in items}
            for future in concurrent.futures.as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:
                    item = futures[future]
                    if ignore_missing and isinstance(exc, dxpy.exceptions.ResourceNotFound):
                        print(f"WARNING: {item} not found, skipping.")
                        continue
                    print(f"Error getting data for {item}: {exc}")
                    raise exc
        return results

    @staticmethod
    def get_report_details(report, project):
        """
        Extracts report details from a given report and project.
        Parameters
        ----------
        report : dict
            dictionary containing report information.
        project : dict
            dictionary containing project information.
        Returns
        -------
        report_details : dict
            dictionary with extracted report details for SNVs and CNV reports.
        """
        report_name = report['describe']['name']
        details = dxpy.bindings.dxdataobject_functions.get_details(report['describe']['id'])
        sample = report_name.split("_")[0]
        assay = project['name'].split("_")[-1]

        report_details = {
            report_name: {
                'project': project['name'],
                'sample': sample,
                'assay': assay,
                'id': report['describe']['id'],
                'clinical_indication': details.get('clinical_indication'),
                'report_type': 'SNV' if 'SNV' in report_name else 'CNV',
                'variants': details.get('included') if 'SNV' in report_name else details.get('variants')
            }
        }
        return report_details

    @staticmethod
    def handle_no_cnv_reports(report_details):
        """
        Handles the case where no CNV reports are found in the report details.
        No CNV reports will return None.

        Parameters
        ----------
        report_details : dict
            dictionary containing report details.
        Returns
        -------
        dict or None
            dictionary of CNV reports or None if no CNV reports found.
        """
        cnv_reports = {
            name: details for name, details in report_details.items()
            if details.get('report_type') == 'CNV'
        }
        if not cnv_reports:
            print("No CNV reports found.")
            return None
        return cnv_reports

    @staticmethod
    def filter_valid_cnv_reports(cnv_reports, excluded_regions_file):
        """
        Filters CNV reports to include only those that are valid based on
        excluded regions. If CNV report has no variants, these reports will
        be checked for excluded regions.

        Parameters
        ----------
        cnv_reports : dict
            Dictionary of CNV reports.
        excluded_regions_file : list
            List of excluded regions files.

        Returns
        -------
        valid_reports : dict
            Dictionary of empty CNV reports.
        """
        valid_reports = {}
        for name, details in cnv_reports.items():
            if details.get('variants', 0) == 0:
                print(f"CNV report {name} has 0 variants.")
                if NMDProcessor.check_excluded_regions_file(excluded_regions_file):
                    valid_reports[name] = details
                else:
                    print(f"Excluded regions was more than 1 for {name}. Skipping.")
        return valid_reports

    @staticmethod
    def check_excluded_regions_file(excluded_regions_file):
        """
        Checks the excluded regions file to determine if it contains
        excluded regions beyond the header.
        Parameters
        ----------
        excluded_regions_file : list
            List of excluded regions files.
        Returns
        -------
        list or None
            Returns the excluded regions file if valid, otherwise None.
        """
        for file in excluded_regions_file:
            file_id = file['id']
            content = dxpy.DXFile(file_id).read().strip().splitlines()
            if len(content) > 1:
                print(f"{file['describe']['name']} contains excluded regions.")
                return None
        print("Excluded regions file has only header. Proceeding.")
        return excluded_regions_file

    @staticmethod
    def filter_overrepresented_indications(report_details, threshold=2):
        """
        Filters out samples with >2 clinical indications in report details.
        Parameters
        ----------
        report_details : dict
            dictionary of report details.
        threshold : int
            Threshold for number of clinical indicatiods
        Returns
        -------
        dict
            Filtered report details with <=2 clinical indications.
        """
        counts = Counter(details['clinical_indication'] for details in report_details.values())
        excluded = {ind for ind, count in counts.items() if count > threshold}
        return {
            name: details for name, details in report_details.items()
            if details['clinical_indication'] not in excluded
        }

# Main processing loop
all_report_details = {}

for proj in projects:
    files = list(dxpy.bindings.search.find_data_objects(classname='file', project=proj['id'], name="^.*NV_\\d+\\.xlsx$", name_mode='regexp', describe=True))
    athena_summary_file = list(dxpy.bindings.search.find_data_objects(classname='file', project=proj['id'], name="^.*1_summary\\.txt$", name_mode='regexp', describe=True))
    excluded_regions_file = list(dxpy.bindings.search.find_data_objects(classname='file', project=proj['id'], name="^.*excluded_intervals.*b38\\.tsv$", name_mode='regexp', describe=True))

    print(f"Found {len(files)} reports in project {proj['name']} ({proj['id']})")

    report_details = NMDProcessor.call_in_parallel(
        func=NMDProcessor.get_report_details,
        items=files,
        project=proj
    )
    report_details = {k: v for report in report_details for k, v in report.items()}
    all_report_details.update(report_details)

    cnv_reports = NMDProcessor.handle_no_cnv_reports(report_details)
    if not cnv_reports:
        continue

    valid_cnv_reports = NMDProcessor.filter_valid_cnv_reports(cnv_reports, excluded_regions_file)