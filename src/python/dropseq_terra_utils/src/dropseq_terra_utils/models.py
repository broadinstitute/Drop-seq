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
from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class WorkflowEmail:
    """
    Represents an email template that should be sent for a specific workflow.
    """
    email_template: str
    workflow_id: str


@dataclass
class FileBytes:
    """
    Represents a file with its name and data.
    """
    filename: str
    data: bytes
    content_type: str = field(default="application/octet-stream")


@dataclass
class WorkspaceInfo:
    """
    Represents a Terra workspace.
    """
    workspace_namespace: str
    workspace_name: str

    @property
    def workspace_url(self) -> str:
        return f"https://app.terra.bio/#workspaces/{self.workspace_namespace}/{self.workspace_name}"


@dataclass
class SubmissionInfo(WorkspaceInfo):
    """
    Represents a Terra submission.
    """
    submission: dict[str, Any]

    @property
    def method_configuration_namespace(self) -> str:
        return self.submission["methodConfigurationNamespace"]

    @property
    def method_configuration_name(self) -> str:
        return self.submission["methodConfigurationName"]

    @property
    def submission_id(self) -> str:
        return self.submission["submissionId"]

    @property
    def submission_url(self) -> str:
        return f"{self.workspace_url}/job_history/{self.submission_id}"

    @property
    def submitter(self) -> str:
        return self.submission["submitter"]

    @property
    def submission_date(self) -> str:
        return self.submission["submissionDate"]

    @property
    def submission_status(self) -> str:
        return self.submission["status"]

    @property
    def is_submission_done(self) -> bool:
        return self.submission_status == "Done"


@dataclass
class WorkflowInfo(SubmissionInfo):
    """
    Represents a Terra workflow.
    """
    submission_workflow: dict[str, Any]
    workflow_metadata: dict[str, Any]

    @property
    def workflow_id(self) -> str:
        return self.submission_workflow["workflowId"]

    @property
    def entity_type(self) -> str:
        return self.submission_workflow["workflowEntity"]["entityType"]

    @property
    def entity_name(self) -> str:
        return self.submission_workflow["workflowEntity"]["entityName"]

    @property
    def workflow_url(self) -> str:
        return f"https://job-manager.dsde-prod.broadinstitute.org/jobs/{self.workflow_id}"


@dataclass
class SubmissionFilters:
    workspaces: list[WorkspaceInfo]
    submission_ids: list[str]
    max_days_old: int


@dataclass
class SmtpSettings:
    """
    Represents the settings needed to send an email via SMTP.
    """
    server: Optional[str]
    port: Optional[int]
    username: Optional[str]
    password: Optional[str]
    tls: bool


@dataclass
class GcloudConfig:
    """
    Represents the configuration needed to interact with Google Cloud.
    """
    name: str
    gcp_project: str
    service_account_email: str
    cloudsdk_config: str

    @property
    def adc_json_path(self) -> str:
        return f"{self.cloudsdk_config}/legacy_credentials/{self.service_account_email}/adc.json"


@dataclass
class HtmlEmailMessage:
    """
    Represents an email message with HTML content.
    """
    email_from: Optional[str] = field(default=None)
    email_to: list[str] = field(default_factory=list)
    subject: Optional[str] = field(default=None)
    body: Optional[str] = field(default=None)
    attachments: list[FileBytes] = field(default_factory=list)
