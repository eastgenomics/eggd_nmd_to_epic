import dxpy
import concurrent.futures
import json
import re
import pandas as pd
import io
from collections import Counter


# Define current project as current workspace
current_project_id = dxpy.WORKSPACE_ID
current_project = dxpy.api.project_describe(current_project_id)
projects = [current_project]


# Class for NMD processing
class NMDProcessor:
    """
    Class for processing NMD reports and filtering based on following criteria:
    - Evaluate if CNV=0 and check excluded regions=0
    - Ignore patients with more than 2 clinical indications
    - Ignore samples where there is no CNV report
    - Extracted reports must have no CNV and SNV variants

    Reports found will be matched to Athena summary files and report outputs written
    in JSON format.
    """
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
                        print(f"{item} not found.")
                        continue
                    print(f"Error for getting data from {item}: {exc}")
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
            Extracted report details for SNVs and CNV reports.
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
        Handles cases where no CNV reports are found in the report details.
        No CNV reports will return None.

        Parameters
        ----------
        report_details : dict
            dictionary containing report details.

        Returns
        -------
        dict or None
            CNV reports or None if no CNV reports found.
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
    def filter_valid_cnv_reports(cnv_reports):
        """
        Filters CNV reports to include only those that are valid based on
        excluded regions. If CNV report has no variants, these reports will
        be checked for excluded regions.

        Parameters
        ----------
        cnv_reports : dict
            Dictionary of CNV reports.

        Returns
        -------
        valid_reports : dict
            CNV reports with no variants and excluded regions.
        """
        valid_reports = {}

        for name, report in cnv_reports.items():
            if report.get("variants", 0) == 0:
                print(f"Checking CNV report: {name}")
                file_id = report.get("id")
                if not file_id:
                    print(f"No file ID found for report: {name}")
                    continue

                try:
                    with dxpy.open_dxfile(file_id, mode='rb') as f:
                        content = f.read()
                        xls = pd.ExcelFile(io.BytesIO(content))
                        if "ExcludedRegions" in xls.sheet_names:
                            df = xls.parse("ExcludedRegions")
                            if not df.empty and len(df.columns) > 0:
                                valid_reports[name] = report
                except Exception as e:
                    print(f"Error processing CNV report {name}: {e}")

        print(f"Valid CNV reports found: {len(valid_reports)}")
        return valid_reports

    @staticmethod
    def filter_overreported_samples(report_details, threshold=2):
        """
        Filters out samples with >2 clinical indications in report details.
        Parameters
        ----------
        report_details : dict
            dictionary of report details.
        threshold : int
            Threshold for number of clinical indications

        Returns
        -------
        dict
            Filtered report details with <=2 clinical indications.
        """
        # Count distinct clinical indications per sample
        sample_to_inds = {}
        for details in report_details.values():
            sample = details.get("sample")
            indication = details.get("clinical_indication")
            if sample is None or indication is None:
                continue
            sample_to_inds.setdefault(sample, set()).add(indication)

        # Identify samples with more than the threshold number of indications
        excluded_samples = {
            sample for sample, inds in sample_to_inds.items()
            if len(inds) > threshold
        }

        # Filter out reports belonging to excluded samples
        return {
            name: details
            for name, details in report_details.items()
            if details.get("sample") not in excluded_samples
        }

    @staticmethod
    def get_reports_with_no_variants(report_details):
        """
        Finds reports with no SNV and CNV variants.
        Parameters
        ----------
        report_details : dict
            dictionary of report details.

        Returns
        -------
        dict
            Reports with no SNV and CNV variants.
        """
        return {
            name: details for name, details in report_details.items()
            if details.get('variants', 0) == 0 and details.get('report_type') in ['SNV', 'CNV']
        }

    @staticmethod
    def get_athena_report(athena_summary_file):
        """
        Finds Athena summary files and gets summary information.
        Parameters
        ----------
        athena_summary_file : list
            List of Athena summary files.
        Returns
        -------
        athena_reports : dict
            Athena report details.
        """
        athena_reports = {}
        for file in athena_summary_file:
            file_id = file['id']
            file_name = file['describe']['name']
            sample_name = file_name.split('_')[0].strip().lower()
            raw_content = dxpy.open_dxfile(file_id, mode='rb').read().strip().splitlines()
            content = [line.decode('utf-8') for line in raw_content]
            athena_reports[sample_name] = {
                "file_id": file_id,
                "content": content
            }
        return athena_reports

    @staticmethod
    def gather_output(filtered_reports, athena_reports):
        """
        Generates output by merging reports and matching Athena summaries.
        Parameters
        ----------
        filtered_reports : dict
            Filtered report details.
        athena_reports : dict
            Athena report details.

        Returns
        -------
        list
            Final output variant reports with matching athena summary
            for NMD cases only.
        """
        final_output = []
        for report_name, details in filtered_reports.items():
            sample = details['sample']
            sample_key = sample.strip().lower()
            athena_summary = athena_reports.get(sample_key)

            output = {
                "report_name": report_name,
                "sample": sample,
                "Epic-InstrumentID" : sample.split('-')[0],
                "Epic-SpecimenID" : sample.split('-')[1],
                "Epic-BatchID" : sample.split('-')[2],
                "project": details['project'],
                "assay": details['assay'],
                "clinical_indication": details['clinical_indication'],
                "report_type": details['report_type'],
                "variants": details['variants'],
                "athena_summary": athena_summary
            }
            final_output.append(output)
        return final_output

# Main processing loop
all_report_details = {}

def main():
    for proj in projects:
        files = list(
            dxpy.bindings.search.find_data_objects(
                classname='file',
                project=proj['id'],
                name="^.*NV_\\d+\\.xlsx$",
                name_mode='regexp',
                describe=True))
        athena_summary_file = list(
            dxpy.bindings.search.find_data_objects(
                classname='file',
                project=proj['id'],
                name="^.*1_summary\\.txt$",
                name_mode='regexp',
                describe=True))

        print(f"Found {len(files)} reports in project {proj['name']} ({proj['id']})")

        report_details = NMDProcessor.call_in_parallel(
            func=NMDProcessor.get_report_details,
            items=files,
            project=proj
        )
        # Get all report details
        report_details = {k: v for report in report_details for k, v in report.items()}
        # Get cases with CNV reports
        cnv_reports = NMDProcessor.handle_no_cnv_reports(report_details)
        if not cnv_reports:
            continue
        # Get CNV reports with no excluded regions
        valid_cnv_reports = NMDProcessor.filter_valid_cnv_reports(cnv_reports)
        # Merge validated CNV reports back with all reports
        merged_reports = {
            name: details for name, details in report_details.items()
            if details.get('report_type') == 'SNV' or name in valid_cnv_reports
        }
        # Get reports with <=2 clinical indications (applies to all report types)
        filtered_reports = NMDProcessor.filter_overreported_samples(merged_reports)
        # Get reports with no CNV and SNV variants
        no_variant_reports = NMDProcessor.get_reports_with_no_variants(filtered_reports)
        # Get Athena summary reports
        athena_reports = NMDProcessor.get_athena_report(athena_summary_file)
        # Create final output
        final_output = NMDProcessor.gather_output(no_variant_reports, athena_reports)
        # Print final output in json format
        for output in final_output:
            print (json.dumps(output, indent = 4))

if __name__ == "__main__":
    main()
