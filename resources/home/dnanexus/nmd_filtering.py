import dxpy
import concurrent.futures
import json
import re
import pandas as pd
import io
import csv
from hl7apy.core import Group
from hl7apy.core import Message
from hl7apy.consts import VALIDATION_LEVEL
from datetime import datetime
from collections import Counter
from utils.other_utils import ReportUtils


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

        # Identify samples with more than threshold and filter out excluded samples
        excluded_samples = {
            sample for sample, inds in sample_to_inds.items()
            if len(inds) > threshold
        }

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

    @staticmethod
    def export_low_coverage(final_output, tsv_path):
        """
        Export reports where the panel coverage at 20x is not 100%.
        Parameters
        ----------
        final_output : list
            List of JSON-like report dicts with athena_summary.
        Returns
        -------
        None
            Writes tsv file with reports having <100% panel coverage.
        """
        rows = []
        for report in final_output:
            athena_summary = report.get("athena_summary", {})
            content = athena_summary.get("content", [])

            for line in content:
                if "of this panel was sequenced to a depth of 20x or greater" in line:
                    # Extract number before % in last line of athena summary
                    line_norm = " ".join(line.split())  # collapse whitespace
                    match = re.search(r"(\d+(?:\.\d+)?)\s*%", line_norm)
                    if match:
                        coverage = float(match.group(1))
                        if coverage < 100:
                            rows.append({
                                "report_name": report.get("report_name", "unknown"),
                                "coverage_percent": coverage
                            })
                    break

        # Write to tsv file
        if rows:
            try:
                with open(tsv_path, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=["report_name", "coverage_percent"], delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows)
            except IOError as e:
                print(f"Error writing to {tsv_path}: {e}")
                raise

            print(f"Exported {len(rows)} reports with <100% coverage to {tsv_path}")
        else:
            print("All reports have 100% coverage.")

    @staticmethod
    def json_to_hl7(json_list):
        """
        Converts a list of JSON report details to HL7 messages.
        Parameters:
        json_list : list
            List of JSON report details.
        Returns:
        list
            List of HL7 messages in ER7 format.
        """
        def pad_obx_to_24_pipes(line):
            """Helper function to make each OBX segment with 24 pipes (25 fields) as expected in EPIC."""
            fields = line.split("|")
            while len(fields) < 25:
                fields.append("")
            return "|".join(fields)

        hl7_messages = []

        for idx, record in enumerate(json_list):
            try:
                msg = Message("ORU_R01", version="2.5.1", validation_level=VALIDATION_LEVEL.TOLERANT)

                specimen_id = record.get("Epic-SpecimenID", "")

                # Handle variant data and create OBX segment
                variants = str(record.get("variants", ""))
                variant_lines = variants.split("\n")

                for i, line in enumerate(variant_lines, start=1):
                    obx_variants = msg.add_segment("OBX")
                    obx_variants.obx_1 = str(i)
                    obx_variants.obx_2 = "ST"
                    obx_variants.obx_3 = f"Variant Genomic details^Variant Genomic details^ATHENA^^^^^^{specimen_id}"
                    obx_variants.obx_5 = line
                    obx_variants.obx_11 = "F"

                # Handle Athena data and create OBX segment
                if record.get("athena_summary"):
                    athena_content = record["athena_summary"].get("content", [])
                    for i, line in enumerate(athena_content, start=1):
                        obx_athena = msg.add_segment("OBX")
                        obx_athena.obx_1 = str(i)
                        obx_athena.obx_2 = "ST"
                        obx_athena.obx_3 = f"Athena Summary^Athena Summary^ATHENA^^^^^^{specimen_id}"
                        obx_athena.obx_5 = line
                        obx_athena.obx_11 = "F"

                # Added metadata in NTE segment
                # ZSP didn't work
                nte = msg.add_segment("NTE")
                nte.nte_3 = f"Project={record.get('project','')}; FileID={record.get('athena_summary',{}).get('file_id','')}"

                # ORU_R01 format will make a header (MSH) automatically
                # so used ER7 format to remove MSH segment
                raw = msg.to_er7()
                processed_lines = []
                for line in raw.split("\r"):
                    if not line.strip() or line.startswith("MSH"):
                        continue
                    if line.startswith("OBX"):
                        line = pad_obx_to_24_pipes(line)
                    processed_lines.append(line)

                raw_no_msh = "\n".join(processed_lines)
                hl7_messages.append(raw_no_msh)

            except Exception as e:
                print(f"Error generating HL7 message: {e}")
                hl7_messages.append(None)

        return hl7_messages

    @staticmethod
    def write_hl7_file(final_output, all_hl7_messages, directory="."):
        """
        Writes all HL7 messages for a run into a single file named:
        <Epic-BatchID>_<project>_NMDs.txt

        Parameters
        ----------
        final_output : list
            List of JSON-like report dicts (each containing Epic-BatchID and project).
        all_hl7_messages : list
            List of HL7 message strings (already generated by json_to_hl7).
        directory : str
            Directory to write the file into (default: current directory).

        Returns
        -------
        str
            Path to the written HL7 file.
        """

        if not final_output:
            raise ValueError("final_output is empty — cannot determine batch/project for filename.")

        # Use the first record to determine run-level metadata
        batch_id = final_output[0].get("Epic-BatchID", "UNKNOWNBATCH")
        project = final_output[0].get("project", "UNKNOWNPROJECT")

        filename = f"{batch_id}_{project}_NMDs.txt"
        filepath = f"{directory}/{filename}"

        try:
            with open(filepath, "w", encoding="utf-8") as f:
                for msg in all_hl7_messages:
                    if msg:
                        f.write(msg + "\n\n")
            return filepath

        except IOError as e:
            print(f"Error writing HL7 file {filepath}: {e}")
            raise

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
        samples_with_valid_cnv = {details['sample'] for details in valid_cnv_reports.values()}
        merged_reports = {
            name: details for name, details in report_details.items()
            if details['sample'] in samples_with_valid_cnv
        }
        # Get reports with <=2 clinical indications (applies to all report types)
        filtered_reports = NMDProcessor.filter_overreported_samples(merged_reports)
        # Get reports with no CNV and SNV variants
        no_variant_reports = NMDProcessor.get_reports_with_no_variants(filtered_reports)
        # Get Athena summary reports
        athena_reports = ReportUtils.get_athena_report(athena_summary_file)
        # Create final output
        final_output = NMDProcessor.gather_output(no_variant_reports, athena_reports)
        # Export poor coverage samples (anything <100% panel coverage)
        NMDProcessor.export_low_coverage(final_output, tsv_path=f"low_coverage_{proj['name']}.tsv")
        # Print final output in json format
        all_hl7_messages = []
        for output in final_output:
            print (json.dumps(output, indent = 4))
            hl7_messages = NMDProcessor.json_to_hl7([output])
            all_hl7_messages.extend(hl7_messages)

            for msg in hl7_messages:
                # Add new line for each segment in hl7 message for readability
                print(msg.replace('\r', '\n'))

        # Make .txt file with all hl7 messages for the run
        filepath = NMDProcessor.write_hl7_file(final_output, all_hl7_messages)
        print(f"Run-level HL7 written to: {filepath}")

        # print how many hl7 messages were made
        print(f"Generated {len(final_output)} HL7 messages for project {proj['name']}.")

        # Print how many were cnvs and how many were snvs
        cnv_count = sum(1 for output in final_output if output['report_type'] == 'CNV')
        snv_count = sum(1 for output in final_output if output['report_type'] == 'SNV')
        print(f"CNV reports: {cnv_count}, SNV reports: {snv_count}")

if __name__ == "__main__":
    main()
