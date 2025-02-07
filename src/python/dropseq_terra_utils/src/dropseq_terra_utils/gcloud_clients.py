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
Services to access Google Cloud related resources including Terra on Google Cloud.
"""
import json
from pathlib import Path
from typing import Any, Optional

import google.auth
import google.cloud.storage
import requests
from google.auth.exceptions import RefreshError
from google.cloud.storage import Blob, Client
from google.oauth2 import service_account
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry as URLLibRetry

try:
    from .models import FileBytes, GcloudConfig
except ImportError:
    from models import FileBytes, GcloudConfig

import google.auth

from google.auth.credentials import TokenState, Credentials
import google.cloud.storage
from google.auth.transport.requests import Request

# Base URL for the FireCloud website
FIRECLOUD_HOST = "https://api.firecloud.org"

# User agent header for requests
USER_AGENT_HEADER = {"User-Agent": "dropseq_terra_utils/0.0.1"}

# Scopes for Google Cloud tokens
GOOGLE_CLOUD_SCOPES = [
    'https://www.googleapis.com/auth/userinfo.profile',
    'https://www.googleapis.com/auth/userinfo.email',
    "https://www.googleapis.com/auth/devstorage.read_only",
]

# Token URI for Google Cloud that isn't present in ADC JSON files, but is needed for the credentials object
GOOGLE_CLOUD_TOKEN_URI = "https://accounts.google.com/o/oauth2/token"


class CredentialsHelper:
    """
    Helper class to manage Google Cloud credentials.
    """

    def __init__(self, gcloud_config: Optional[GcloudConfig] = None):
        self.credentials: Credentials

        if gcloud_config:
            with open(gcloud_config.adc_json_path) as file:
                credentials_info = json.load(file)
            if not isinstance(credentials_info, dict):
                raise ValueError(f"ADC JSON file must contain a JSON object: {gcloud_config.adc_json_path}")
            if "token_uri" not in credentials_info:
                credentials_info["token_uri"] = GOOGLE_CLOUD_TOKEN_URI
            self.credentials = service_account.Credentials.from_service_account_info(credentials_info,
                                                                                     scopes=GOOGLE_CLOUD_SCOPES)
        else:
            self.credentials, _ = google.auth.default(scopes=GOOGLE_CLOUD_SCOPES)

        self.project: Optional[str] = gcloud_config.gcp_project if gcloud_config else None
        self.request: Optional[Request] = None

    def fresh_credentials(self) -> Credentials:
        """
        Refresh the credentials if they are not fresh.
        """
        if self.credentials.token_state != TokenState.FRESH:
            if not self.request:
                self.request = Request()
            self.credentials.refresh(self.request)
        return self.credentials

    def is_valid_credentials(self) -> bool:
        """
        Check if the credentials are valid quietly, without raising an exception.
        """
        try:
            self.fresh_credentials()
        except RefreshError:
            return False
        return True

    def make_auth_header(self) -> dict[str, str]:
        """
        Create an authorization header for a request with a fresh token.
        """
        return {"Authorization": f"Bearer {self.fresh_credentials().token}"}


class GcsClient:
    """
    Client to access Google Cloud Storage.
    """

    def __init__(self, credentials: CredentialsHelper, max_file_size: int = -1):
        self.credentials: CredentialsHelper = credentials

        if self.credentials.project:
            self.client: Client = google.cloud.storage.Client(
                credentials=self.credentials.fresh_credentials(),
                project=self.credentials.project,
            )
        else:
            # Use the project found by the environment
            self.client: Client = google.cloud.storage.Client(
                credentials=self.credentials.fresh_credentials(),
            )
        self.max_file_size: int = max_file_size

    def download_file(self, gcs_path: str) -> FileBytes:
        """
        Download a file from Google Cloud Storage.
        """
        blob = Blob.from_string(gcs_path, self.client)
        blob.reload()
        size = blob.size
        if self.max_file_size > -1:
            if self.max_file_size < size:
                raise ValueError(f"File size is too large to download: {size} bytes")
        return FileBytes(Path(blob.name).name, blob.download_as_bytes(), blob.content_type)


class TerraClient:
    """
    Client to access Terra on Google Cloud.
    """

    def __init__(self, credentials: CredentialsHelper):
        self.firecloud_host: str = FIRECLOUD_HOST
        self.credentials: CredentialsHelper = credentials

    def _make_request(self, url: str) -> Any:
        """
        Make a request to the Terra API with retries for common errors.

        :param url: The URL to make the request to.
        :return: The content parsed from json.
        """
        retry_strategy = URLLibRetry(
            total=3,
            backoff_factor=1.2,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["HEAD", "GET", "OPTIONS", "POST"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount("https://", adapter)
        headers = {}
        headers.update(self.credentials.make_auth_header())
        headers.update(USER_AGENT_HEADER)
        response = session.get(url, headers=headers)
        return json.loads(response.text)

    def get_user_details(self) -> dict[str, Any]:
        """
        Get the current Terra user.

        :return: The current Terra user parsed from json.
        """
        url = f"{self.firecloud_host}/me?userDetailsOnly=true"
        return self._make_request(url)

    def get_workspaces(self) -> list[dict[str, Any]]:
        """
        Get a list of workspaces from Terra.

        :return: A list of workspaces parsed from json.
        """
        url = f"{self.firecloud_host}/api/workspaces"
        return self._make_request(url)

    def get_submissions(self, workspace_namespace: str, workspace_name: str) -> list[dict[str, Any]]:
        """
        Get a list of submissions for a workspace, without submission workflow information.

        :param workspace_namespace: The namespace of the workspace.
        :param workspace_name: The name of the workspace.
        :return: A list of submissions parsed from json.
        """
        url = f"{self.firecloud_host}/api/workspaces/{workspace_namespace}/{workspace_name}/submissions"
        return self._make_request(url)

    def get_submission(self, workspace_namespace: str, workspace_name: str, submission_id: str) -> dict[str, Any]:
        """
        Get a submission information, including submission workflows.

        :param workspace_namespace: The namespace of the workspace.
        :param workspace_name: The name of the workspace.
        :param submission_id: The id of the submission.
        :return: A submission parsed from json.
        """
        url = f"{self.firecloud_host}/api/workspaces/{workspace_namespace}/{workspace_name}/submissions/{submission_id}"
        return self._make_request(url)

    def get_workflow(self, workflow_id: str) -> dict[str, Any]:
        """
        Get workflow metadata, including subworkflow metadata.

        :param workflow_id: The id of the workflow.
        :return: Workflow metadata parsed from json.
        """
        url = f"{self.firecloud_host}/api/workflows/v1/{workflow_id}/metadata?expandSubWorkflows=true"
        return self._make_request(url)
