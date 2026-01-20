import dxpy

class ReportUtils:
    """
    Class for gathering sample metadata for PASS variants
    and NMD samples
    """
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
            with dxpy.open_dxfile(file_id, mode='rb') as f:
                raw_content = f.read().strip().splitlines()
                content = [line.decode('utf-8') for line in raw_content]
            athena_reports[sample_name] = {
                "file_id": file_id,
                "file_name": file_name,
                "content": content
            }
        return athena_reports

    def pad_obx_to_24_pipes(line):
        """Helper function to make each OBX segment with 24 pipes (25 fields) as expected in EPIC."""
        fields = line.split("|")
        while len(fields) < 25:
            fields.append("")
        return "|".join(fields)