#!/usr/bin/env python3
# MIT License
#
# Copyright 2025 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Email templates for sending emails based on workflow outputs.
"""
from abc import abstractmethod, ABC
from io import StringIO
from typing import Any, Optional, Union

import pandas as pd
from pandas import DataFrame

try:
    from .models import FileBytes, HtmlEmailMessage, WorkflowInfo
except ImportError:
    from models import FileBytes, HtmlEmailMessage, WorkflowInfo


class EmailTemplate(ABC):
    """
    An email template for sending emails based on workflow outputs.
    """

    def __init__(self, config_name: str):
        self.config_name: str = config_name

    @property
    @abstractmethod
    def workflow_names(self) -> list[str]:
        """
        Get the names of the workflows that this email template matches.
        This list should include all prior names of the workflow that have been used.
        """
        pass

    @abstractmethod
    def find_output_paths(self, workflow_info: WorkflowInfo) -> dict[str, Optional[str]]:
        """
        Find the output paths for the email template within the workflow data.

        The keys of the dictionary will be the same as those passed to the make_message method.
        The values should be the paths to the output files within the workflow data, or None if the output is not found.
        """
        pass

    @abstractmethod
    def make_message(self, workflow_info: WorkflowInfo, output_data: dict[str, FileBytes]) -> HtmlEmailMessage:
        """
        Create the email message based on the workflow data and output data.
        """
        pass

    @staticmethod
    def _tsv_to_df(tsv_bytes: FileBytes) -> DataFrame:
        """
        Convert a TSV file to a DataFrame.
        """
        return pd.read_csv(StringIO(tsv_bytes.data.decode("utf-8")), sep="\t", dtype=str)

    @staticmethod
    def _find_call_output(
            workflow_info: WorkflowInfo,
            call_name: Union[str, list[str]],
            output_name: Union[str, list[str]],
            skip_call_cache: bool = False
    ) -> Optional[Any]:
        """
        Find the output for a call.
        """
        if not isinstance(call_name, list):
            call_name = [call_name]
        if not isinstance(output_name, list):
            output_name = [output_name]
        call_success = EmailTemplate.__find_call_success(workflow_info.workflow_metadata, call_name, skip_call_cache)
        return EmailTemplate.__get_call_output(call_success, output_name)

    @staticmethod
    def _subject_line(message_type: str, workflow_info: WorkflowInfo) -> str:
        return f"{message_type} for {workflow_info.entity_name}"

    @staticmethod
    def _body_header(message_type: str, workflow_info: WorkflowInfo) -> str:
        return f"""
        <H2>
            {message_type} for
        </H2>
        <table border="0">
            <tr>
                <td>
                    <b>{workflow_info.entity_type}</b>
                </td>
                <td>
                    <a href="{workflow_info.workflow_url}">
                        {workflow_info.entity_name}
                    </a>
                </td>
            </tr>
            <tr>
                <td>
                    <b>submission</b>
                </td>
                <td>
                    <a href="{workflow_info.submission_url}">
                        {workflow_info.method_configuration_namespace}/{workflow_info.method_configuration_name}
                    </a>
                </td>
            </tr>
            <tr>
                <td>
                    <b>workspace</b>
                </td>
                <td>
                    <a href="{workflow_info.workspace_url}">
                        {workflow_info.workspace_namespace}/{workflow_info.workspace_name}
                    </a>
                </td>
            </tr>
            <tr>
                <td>
                    <b>date submitted</b>
                </td>
                <td>
                    {workflow_info.submission_date}
                </td>
            </tr>
        </table>
        """

    @staticmethod
    def __find_call_success(
            workflow_metadata: dict[str, Any],
            call_names: list[str],
            skip_call_cache: bool = False
    ) -> dict[str, Any]:
        """
        Get the call that matches the passed in name.
        """
        for call_name in call_names:
            call = workflow_metadata.get("calls", {}).get(call_name)
            if call:
                for attempt in call:
                    # Look for the first entry with a return code of zero
                    if attempt.get("returnCode") == 0:
                        if skip_call_cache and attempt.get("callCaching", {}).get("hit"):
                            return {}
                        else:
                            return attempt
        return EmailTemplate.__find_subworkflow_call_success(workflow_metadata, call_names, skip_call_cache)

    @staticmethod
    def __find_subworkflow_call_success(
            workflow_metadata: dict[str, Any],
            call_names: list[str],
            skip_call_cache: bool = False
    ) -> dict:
        """
        Find the successful call by looking in subworkflow metadata.
        """
        subworkflows = workflow_metadata.get("calls", {})
        for subworkflow_name in subworkflows:
            subworkflow = subworkflows[subworkflow_name]
            for attempt in subworkflow:
                if "subWorkflowMetadata" not in attempt:
                    continue
                call_success = EmailTemplate.__find_call_success(
                    attempt["subWorkflowMetadata"],
                    call_names,
                    skip_call_cache
                )
                if call_success:
                    return call_success
        return {}

    @staticmethod
    def __get_call_output(call_success: dict[str, Any], output_names: list[str]) -> Optional[Any]:
        """
        Get the output for a call success.
        """
        if "outputs" not in call_success:
            return None
        outputs = call_success["outputs"]
        for output_name in output_names:
            if output_name in call_success["outputs"]:
                return outputs[output_name]
        return None


class PdfTearsheet(EmailTemplate, ABC):
    """
    An email template for sending tearsheets with a PDF attachment.
    """

    _TEARSHEET_PDF_KEY = "tearsheet_pdf"

    def __init__(self, config_name: str, message_type: str):
        super().__init__(config_name)
        self.message_type: str = message_type

    def make_message(self, workflow_info: WorkflowInfo, output_data: dict[str, FileBytes]) -> HtmlEmailMessage:
        email_message = HtmlEmailMessage()
        email_message.subject = self._subject_line(self.message_type, workflow_info)
        email_message.body = self._body_html(workflow_info, output_data)
        email_message.attachments = [output_data[self._TEARSHEET_PDF_KEY]]
        return email_message

    def _body_html(self, workflow_info: WorkflowInfo, output_data: dict[str, FileBytes]) -> str:
        return f"""\n<div>\n{self._body_header(self.message_type, workflow_info)}</div>\n"""


class TsvTearsheet(PdfTearsheet, ABC):
    """
    An email template for sending tearsheets with PDF attachments and an email body generated from a TSV.
    """

    _TEARSHEET_TSV_KEY = "tearsheet_tsv"

    def __init__(self, config_name: str, message_type: str):
        super().__init__(config_name, message_type)

    def _body_html(self, workflow_info: WorkflowInfo, output_data: dict[str, FileBytes]) -> str:
        tearsheet_df = EmailTemplate._tsv_to_df(output_data[self._TEARSHEET_TSV_KEY])
        body = f"""<div>\n{self._body_header(self.message_type, workflow_info)}\n<table border="1">\n"""
        for _, row in tearsheet_df.iterrows():
            body += f"""<tr><td>{row['label']}</td><td>{row['value']}</td></tr>\n"""
        body += """</table>\n</div>\n"""
        return body


class CbrbTearsheet(TsvTearsheet):
    def __init__(self):
        super().__init__("cbrb_tearsheet", "CBRB report")

    @property
    def workflow_names(self) -> list[str]:
        return [
            "dropseq_cbrb",
            "optimus_cbrb",
            "optimus_dropulation",
            "optimus_auto_dropulation",
            "optimus_auto_cell_selection",
        ]

    def find_output_paths(self, workflow_info: WorkflowInfo) -> dict[str, Optional[str]]:
        return {
            self._TEARSHEET_TSV_KEY:
                EmailTemplate._find_call_output(
                    workflow_info,
                    "dropseq_cbrb.make_cbrb_0_3_0_tear_sheet_properties",
                    "out_file",
                ),
            self._TEARSHEET_PDF_KEY:
                EmailTemplate._find_call_output(
                    workflow_info,
                    "dropseq_cbrb.cbrb_0_3_0_tear_sheet",
                    "out_file",
                ),
        }


class CellSelectionTearsheet(TsvTearsheet):
    """
    An email template for sending Cell Selection tearsheets.
    """

    def __init__(self):
        super().__init__("cell_selection_tearsheet", "CallSTAMPs report")

    @property
    def workflow_names(self) -> list[str]:
        return [
            "cell_selection",
            "selection_dropulation",
            "optimus_dropulation",
            "optimus_auto_dropulation",
            "optimus_auto_cell_selection",
            "manual_cell_selection_dropulation",
        ]

    def find_output_paths(self, workflow_info: WorkflowInfo) -> dict[str, Optional[str]]:
        return {
            self._TEARSHEET_TSV_KEY:
                EmailTemplate._find_call_output(
                    workflow_info,
                    "cell_selection.make_standard_analysis_tear_sheet",
                    "out_file",
                ),
            self._TEARSHEET_PDF_KEY:
                EmailTemplate._find_call_output(
                    workflow_info,
                    "cell_selection.call_stamps",
                    "out_pdf",
                ),
        }


class DropulationTearsheet(PdfTearsheet):
    """
    An email template for sending Dropulation tearsheets.
    """

    def __init__(self):
        super().__init__("dropulation_tearsheet", "Dropulation Tear Sheet")

    @property
    def workflow_names(self) -> list[str]:
        return [
            "dropulation",
            "standard_analysis",
            "selection_dropulation",
            "optimus_dropulation",
            "optimus_auto_dropulation",
            "manual_cell_selection_dropulation",
        ]

    def find_output_paths(self, workflow_info: WorkflowInfo) -> dict[str, Optional[str]]:
        return {
            self._TEARSHEET_PDF_KEY:
                EmailTemplate._find_call_output(
                    workflow_info,
                    "dropulation.donor_assignment_qc",
                    "out_tear_sheet_pdf",
                ),
        }


def email_templates_list() -> list[EmailTemplate]:
    """
    Get a list of all available email templates.
    """
    return [
        CbrbTearsheet(),
        CellSelectionTearsheet(),
        DropulationTearsheet(),
    ]
