import dxpy
from hl7apy.core import Message
from hl7apy.consts import VALIDATION_LEVEL
from datetime import datetime
from utils.other_utils import ReportUtils


class PassProcessor:
    """
    Class for processing PASS SNV reports and generating HL7 messages with Athena summaries.
    """

    @staticmethod
    def get_filtered_reports(snv_reports_files):
        """
        Build filtered_reports dict with metadata from SNV report files.
        Parameters
        ----------
        snv_reports_files : list
            List of SNV report files.
        Returns
        -------
        filtered reports : dict
            Filtered reports with metadata.
        """
        filtered_reports = {}
        for file in snv_reports_files:
            desc = file['describe']
            file_name = desc.get('name', '')
            sample = file_name.split('_')[0].strip()
            project_id = desc.get('project', '')

            # Resolve project name from ID
            project_name = ''
            if project_id:
                try:
                    project_name = dxpy.api.project_describe(project_id).get('name', '')
                except Exception as e:
                    print(f"Warning: could not resolve project name for {project_id}: {e}")

            filtered_reports[file_name] = {
                "sample": sample,
                "project_id": project_id,
                "project_name": project_name,
                "assay": desc.get('assay', ''),
                "clinical_indication": desc.get('clinical_indication', ''),
                "report_type": desc.get('report_type', ''),
                "variants": desc.get('variants', [])
            }
        return filtered_reports

    @staticmethod
    def gather_all_outputs(filtered_reports, athena_reports):
        """
        Merge SNV reports with matching Athena summaries.
        Parameters
        ----------
        filtered_reports : dict
            Filtered SNV reports.
        athena_reports : dict
            Athena summaries keyed by sample.

        Returns
        -------
        final_output : list
            Merged athena output list.
        """
        final_output = []

        for report_name, details in filtered_reports.items():
            sample = details['sample']
            sample_key = sample.strip().lower()

            # Find athena files (case insensitive)
            athena_summary = None
            for k, v in athena_reports.items():
                if k.strip().lower() == sample_key:
                    athena_summary = v
                    break

            # Print result to make sure samples matched are same
            if athena_summary:
                print(f"Matched sample for SNV: {sample}")
                print(f"Athena report matched: file_name={athena_summary.get('file_name')}")
            else:
                print(f"WARNING: No Athena summary file found for sample '{sample}' (key: {sample_key})")

            parts = sample.split('-') if sample else []
            epic_instrument_id = parts[0] if len(parts) > 0 else ""
            epic_specimen_id = parts[1] if len(parts) > 1 else ""
            epic_batch_id = parts[2] if len(parts) > 2 else ""

            output = {
                "report_name": report_name,
                "sample": sample,
                "Epic-InstrumentID": epic_instrument_id,
                "Epic-SpecimenID": epic_specimen_id,
                "Epic-BatchID": epic_batch_id,
                "project_id": details.get('project_id', ""),
                "project_name": details.get('project_name', ""),
                "assay": details.get('assay', ""),
                "clinical_indication": details.get('clinical_indication', ""),
                "report_type": details.get('report_type', ""),
                "variants": details.get('variants', []),
                "athena_summary": athena_summary
            }
            final_output.append(output)

        return final_output

    @staticmethod
    def generate_pass_hl7_message(file, athena_summary, project_name):
        """
        Generate HL7 message for SNV file and Athena summary.
        Parameters
        ----------
        file : dict
            SNV report file's metadata.
        athena_summary : dict
            Athena summary data.
        project_name : str
            Project name.
        Returns
        -------
        hl7_message : str
            HL7 message with data interpretation location and athena summary.
        """
        try:
            msg = Message("ORU_R01", version="2.5.1", validation_level=VALIDATION_LEVEL.TOLERANT)

            specimen_id = file.get("sample", "")

            # Create NTE with directory path based on project name
            base_path = ""
            # Using '\\' to avoid python escape issues and for HL7 compatibility
            if project_name and "CEN" in project_name.upper():
                base_path = "\\\\clingen\\cg\\Regional Genetics Laboratories\\Molecular Genetics\\Data archive\\Sequencing HT\\CEN\\Run folders"
            elif project_name and ("WES" in project_name.upper() or "TWE" in project_name.upper()):
                base_path = "\\\\clingen\\cg\\Regional Genetics Laboratories\\Molecular Genetics\\Data archive\\Sequencing HT\\WES"

            cleaned_project_name = project_name.replace("002_", "") if project_name else ""
            directory_path = f"{base_path}\\{cleaned_project_name}"

            # OBX segments for Athena summary
            content_lines = athena_summary.get("content", []) if athena_summary else []
            content_lines = [line.strip() for line in content_lines if line.strip()]

            for i, line in enumerate(content_lines, start=1):
                obx = msg.add_segment("OBX")
                obx.obx_1 = str(i)
                obx.obx_2 = "TX"
                obx.obx_3 = f"Athena Summary^Athena Summary^ATHENA^^^^^^{specimen_id}"
                obx.obx_5 = line
                obx.obx_11 = "F"
            # Add NTE segment after OBX segments
            nte = msg.add_segment("NTE")
            nte.nte_3 = directory_path

            raw = msg.to_er7()

            processed_lines = []
            for line in raw.split("\r"):
                if not line.strip() or line.startswith("MSH"):
                    continue
                if line.startswith("OBX"):
                    line = ReportUtils.pad_obx_to_24_pipes(line)
                processed_lines.append(line)

            return "\n".join(processed_lines)

        except Exception as e:
            print(f"Error generating HL7 message: {e}")
            return None

    @staticmethod
    def write_pass_hl7_file(final_output, all_hl7_messages, directory="."):
        """
        Writes all PASS HL7 messages into:
        <Epic-BatchID>_<project>_pass_samples.txt
        """

        if not final_output:
            raise ValueError("final_output is empty.")

        batch_id = final_output[0].get("Epic-BatchID", "UNKNOWNBATCH")
        project = final_output[0].get("project_name", "UNKNOWNPROJECT")

        filename = f"{batch_id}_{project}_pass_samples.txt"
        filepath = f"{directory}/{filename}"

        try:
            with open(filepath, "w", encoding="utf-8") as f:
                for msg in all_hl7_messages:
                    if msg:
                        f.write(msg + "\n\n")
            return filepath

        except IOError as e:
            print(f"Error writing PASS HL7 file {filepath}: {e}")
            raise

def main():
    # Define current project as current workspace
    current_project_id = dxpy.WORKSPACE_ID
    current_project = dxpy.api.project_describe(current_project_id)
    projects = [current_project]

    for proj in projects:
        snv_reports_files = list(
            dxpy.bindings.search.find_data_objects(
                classname='file',
                name=".*SNV_.*\\.xlsx$",
                name_mode='regexp',
                describe=True,
                project=proj['id']
            )
        )
        athena_summary_file = list(
            dxpy.bindings.search.find_data_objects(
                classname='file',
                project=proj['id'],
                name="^.*1_summary\\.txt$",
                name_mode='regexp',
                describe=True
            )
        )

    print("Found", len(snv_reports_files), "SNV report files")

    all_filenames = [file['describe']['name'] for file in snv_reports_files]
    print("SNV report files found:", all_filenames)

    # Build filtered_reports from SNV files
    filtered_reports = PassProcessor.get_filtered_reports(snv_reports_files)
    # Get Athena summaries
    athena_reports = ReportUtils.get_athena_report(athena_summary_file)
    # Merge into final output
    final_pass_output = PassProcessor.gather_all_outputs(filtered_reports, athena_reports)
    # Generate HL7 message for each pass sample
    all_hl7_messages = []
    hl7_count = 0
    for output in final_pass_output:
        filename = output["report_name"]
        print("Processing:", filename)

        hl7_message = PassProcessor.generate_pass_hl7_message(
            {"describe": {"name": filename}, "sample": output["sample"]},
            output["athena_summary"],
            output["project_name"]
        )

        if hl7_message:
            # Fix HL7 escape sequences for backslashes (\ is escape character in hl7)
            hl7_message = hl7_message.replace("\\E\\", "\\")
            print("HL7 message generated for", filename)
            print(hl7_message.replace('\r', '\n'))
            # Add to all messages list
            all_hl7_messages.append(hl7_message)
            hl7_count += 1
        else:
            print("No HL7 message for", filename)

    # Write all HL7 messages to file
    filepath = PassProcessor.write_pass_hl7_file(final_pass_output, all_hl7_messages)
    print(f"PASS HL7 file written to: {filepath}")

    # Print how many hl7 messages were created from count
    print(f"\nTotal HL7 messages generated: {hl7_count}")

if __name__ == "__main__":
    main()
