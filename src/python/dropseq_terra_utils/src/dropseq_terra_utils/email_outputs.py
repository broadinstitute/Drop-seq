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
Email outputs from Terra workflows.
"""
import argparse
import html
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

import yaml
from dateutil import parser as date_parser

try:
    from . import cli
    from . import data_store
    from . import email_templates
    from .data_store import DataStore
    from .email_clients import EmailClient, SmtpEmailClient
    from .gcloud_clients import CredentialsHelper, GcsClient, TerraClient, GcloudClients
    from .email_templates import EmailTemplate
    from .models import WorkflowEmail, SubmissionInfo, WorkflowInfo, HtmlEmailMessage, SmtpSettings, WorkspaceInfo, \
        SubmissionFilters, GcloudConfig
except ImportError:
    import cli
    import data_store
    import email_templates
    from data_store import DataStore
    from email_clients import EmailClient, SmtpEmailClient
    from gcloud_clients import CredentialsHelper, GcsClient, TerraClient, GcloudClients
    from email_templates import EmailTemplate
    from models import WorkflowEmail, SubmissionInfo, WorkflowInfo, HtmlEmailMessage, SmtpSettings, WorkspaceInfo, \
        SubmissionFilters, GcloudConfig


class SystemClients:
    """
    Collection of clients for various system services.
    """

    def __init__(self, emailer: EmailClient, ds: DataStore):
        self.email_client: EmailClient = emailer
        self.data_store: DataStore = ds


def add_subparser(subparsers) -> None:
    parser = subparsers.add_parser("email_outputs", description=__doc__)
    add_arguments(parser)


def add_arguments(parser) -> None:
    parser.add_argument(
        "--project-metadata", "-p",
        help="Project metadata YAML with service accounts to use.  Default: use the default service account.",
    )
    parser.add_argument(
        "--workspace", "-w", action="append",
        help="Workspace namespace/name(s) to process.  Default: all non-public workspaces",
    )
    parser.add_argument(
        "--submission", "-s", action="append",
        help="Submission ID(s) to process.  Default: all submissions in the workspace(s)",
    )
    parser.add_argument(
        "--max-days-old", type=int, default=21,
        help="Maximum number of days old a submission can be to process.  Default: %(default)s",
    )
    parser.add_argument(
        "--max-gcs-bytes", type=int, default=1_000_000_000,
        help="Maximum file size in bytes to download from GCS.  Default: %(default)s",
    )
    parser.add_argument(
        "--data-store", "-d", default=data_store.default_data_store_dir(),
        help="Directory to store files such as tracking prior progress.  Default: %(default)s",
    )
    parser.add_argument(
        "--recheck", "-r", type=int, default=300,
        help="Number of seconds to wait between rechecking for email outputs.  Default: %(default)s",
    )
    parser.add_argument("--smtp-username", help="Optional SMTP username.")
    parser.add_argument("--smtp-password", help="Optional SMTP password.")
    parser.add_argument("--smtp-server", default="smtp", help="SMTP server.  Default: %(default)s")
    parser.add_argument("--smtp-port", type=int, default=25, help="SMTP port.  Default: %(default)s")
    parser.add_argument("--smtp-tls", action="store_true", help="Use TLS for SMTP. Default: %(default)s")
    parser.add_argument("--email-from", help="Default email from address.  Default: derived from the user and host")
    parser.add_argument("--errors-to", help="Email address to send errors to.")


def run(options):
    cli.logger.info("Starting")

    if options.project_metadata:
        gcloud_configs = get_service_accounts(options.project_metadata)
        credential_helpers = get_terra_credentials(gcloud_configs)
    else:
        cli.logger.info("Using application default credentials")
        credential_helpers = [CredentialsHelper()]

    smtp_settings = SmtpSettings(
        options.smtp_server,
        options.smtp_port,
        options.smtp_username,
        options.smtp_password,
        options.smtp_tls,
    )
    smtp = SmtpEmailClient(smtp_settings, options.email_from, options.errors_to)
    ds = DataStore(options.data_store)
    system_clients = SystemClients(smtp, ds)

    gcloud_clients_list = []
    for credentials_helper in credential_helpers:
        gcs = GcsClient(credentials_helper, options.max_gcs_bytes)
        terra = TerraClient(credentials_helper)
        gcloud_clients_list.append(GcloudClients(gcs, terra))

    email_templates_list = email_templates.email_templates_list()

    submission_filters = SubmissionFilters(
        [WorkspaceInfo(*workspace.split("/", 1)) for workspace in (options.workspace or [])],
        options.submission or [],
        options.max_days_old,
    )

    run_with_retries(
        lambda: send_workspace_emails(system_clients, gcloud_clients_list, email_templates_list, submission_filters),
        system_clients,
        options.recheck,
    )

    cli.logger.info("Done")
    return 0


def get_service_accounts(project_metadata_path: str) -> list[GcloudConfig]:
    """
    Get the service accounts from a project metadata YAML.

    :param project_metadata_path: The project YAML file path.
    :return: The service accounts.
    """
    with open(project_metadata_path, 'r') as file:
        yaml_contents = yaml.safe_load(file)
    gcloud_configs = []
    for project in yaml_contents["projects"]:
        name = project["name"]
        gcp_project = project.get("gcp_project", name)
        service_account_email = project["service_account"]
        cloudsdk_config = project["cloudsdk_config"]
        gcloud_configs.append(GcloudConfig(name, gcp_project, service_account_email, cloudsdk_config))
    return gcloud_configs


def get_terra_credentials(gcloud_configs: list[GcloudConfig]) -> list[CredentialsHelper]:
    """
    Retrieve the credentials for services accounts that are also Terra users.
    """
    credentials_helpers = []
    emails: dict[str, GcloudConfig] = {}

    for gcloud_config in gcloud_configs:
        email = gcloud_config.service_account_email
        if email in emails:
            cli.logger.info(f"Terra user {email} using config {emails[email].name} instead of {gcloud_config.name}")
            continue
        if not Path(gcloud_config.adc_json_path).exists():
            cli.logger.info(f"ADC JSON file for {email} not found via {gcloud_config.name}")
            continue

        credentials_helper = CredentialsHelper(gcloud_config)
        if not credentials_helper.is_valid_credentials():
            cli.logger.info(f"Failed to load credentials for {email} via {gcloud_config.name}")
            continue

        terra_client = TerraClient(credentials_helper)
        user_details = terra_client.get_user_details()
        if not isinstance(user_details, dict) or not user_details.get("enabled", False):
            cli.logger.info(f"Terra user {email} found not for project {gcloud_config.name}")
            continue

        cli.logger.info(f"Using credentials {email} via project {gcloud_config.name}")
        credentials_helpers.append(credentials_helper)
        emails[email] = gcloud_config

    return credentials_helpers


def get_workspaces(gcloud_clients: GcloudClients, submission_filters: SubmissionFilters) -> list[WorkspaceInfo]:
    """
    Get the workspaces to process.

    :param gcloud_clients: The Google Cloud clients.
    :param submission_filters: The workflow filters.
    :return: The list of workspaces to process.
    """
    if submission_filters.workspaces:
        workspace_infos = submission_filters.workspaces
    else:
        all_workspaces = gcloud_clients.terra.get_workspaces()
        filtered_workspaces = [
            workspace["workspace"]
            for workspace in all_workspaces
            if (not workspace["public"]
                and workspace["accessLevel"] != "NO ACCESS"
                and workspace["workspace"]["cloudPlatform"] == "Gcp")
        ]
        sorted_workspaces = sorted(filtered_workspaces, key=lambda w: w["lastModified"], reverse=True)
        workspace_infos = [
            WorkspaceInfo(workspace["namespace"], workspace["name"])
            for workspace in sorted_workspaces
        ]
    return workspace_infos


def get_submissions(
        gcloud_clients: GcloudClients,
        workspace: WorkspaceInfo,
        submission_filters: SubmissionFilters,
) -> list[SubmissionInfo]:
    """
    Get the submissions to process.

    :param gcloud_clients: The Google Cloud clients.
    :param workspace: The workspace to process.
    :param submission_filters: The workflow filters.
    :return: The list of submissions to process.
    """
    all_submissions = gcloud_clients.terra.get_submissions(workspace.workspace_namespace, workspace.workspace_name)

    submission_infos = [
        SubmissionInfo(workspace.workspace_namespace, workspace.workspace_name, submission)
        for submission in sorted(
            all_submissions,
            key=lambda s: s["submissionDate"],
        )
    ]

    if submission_filters.submission_ids:
        submission_infos = [
            submission
            for submission in submission_infos
            if submission.submission_id in submission_filters.submission_ids
        ]

    if submission_filters.max_days_old > -1:
        submission_infos = [
            submission
            for submission in submission_infos
            if (
                       datetime.now(timezone.utc) - date_parser.isoparse(submission.submission_date)
               ).days <= submission_filters.max_days_old
        ]

    return submission_infos


def get_submission_workflow_names(submission_detail: dict[str, Any]) -> list[str]:
    """
    Get the probable workflow names from a submission detail.

    Whenever a method configuration is changed the prior configuration is both marked as deleted and renamed with a new
    suffix.

    So we can't always retrieve the method configuration from the API to find out what workflow was used.
    That exact information is also available in the workflow metadata, but we don't have that yet.
    Instead, we can look at the input resolutions to see what workflow names were used.
    Most likely, the workflow name is the first part of the input names. Some call names may slip in.

    :param submission_detail: The submission detail.
    :return: The list of probable workflow names.
    """
    # Get the input resolutions from the submission detail
    input_resolutions: list[dict[str, Any]] = submission_detail["workflows"][0]["inputResolutions"]
    # Split on the first "." in the input name to get the workflow name, or sometimes the call name
    workflow_names: list[str] = [res["inputName"].split(".", 1)[0] for res in input_resolutions]
    # Count the occurrences of each workflow name
    workflow_counts: dict[str, int] = {name: workflow_names.count(name) for name in set(workflow_names)}
    # Sort the workflow names by the count of occurrences
    return sorted(workflow_counts, key=workflow_counts.get, reverse=True)


def send_workspace_emails(
        system_clients: SystemClients,
        gcloud_clients_list: list[GcloudClients],
        email_templates_list: list[EmailTemplate],
        submission_filters: SubmissionFilters,
) -> None:
    for gcloud_clients in gcloud_clients_list:
        workspaces = get_workspaces(gcloud_clients, submission_filters)
        for workspace in workspaces:
            process_workspace(system_clients, gcloud_clients, email_templates_list, workspace, submission_filters)


def process_workspace(
        system_clients: SystemClients,
        gcloud_clients: GcloudClients,
        email_templates_list: list[EmailTemplate],
        workspace: WorkspaceInfo,
        submission_filters: SubmissionFilters,
) -> None:
    cli.logger.debug(f"Processing workspace {workspace.workspace_namespace}/{workspace.workspace_name}")
    submission_infos = get_submissions(gcloud_clients, workspace, submission_filters)
    for submission_info in submission_infos:
        if system_clients.data_store.has_done_submission(submission_info.submission_id):
            cli.logger.debug(f"Submission {submission_info.submission_id} was already done")
            continue

        cli.logger.debug(f"Processing submission {submission_info.submission_id}")
        process_submission(system_clients, gcloud_clients, email_templates_list, submission_info)

        # If the submission is done, mark it as done and don't process it again
        if submission_info.is_submission_done:
            cli.logger.info(f"Marking submission {submission_info.submission_id} as done")
            system_clients.data_store.add_done_submission(submission_info.submission_id)


def process_submission(
        system_clients: SystemClients,
        gcloud_clients: GcloudClients,
        email_templates_list: list[EmailTemplate],
        submission_info: SubmissionInfo,
) -> None:
    submission_detail = gcloud_clients.terra.get_submission(
        submission_info.workspace_namespace,
        submission_info.workspace_name,
        submission_info.submission_id,
    )
    maybe_workflow_names = get_submission_workflow_names(submission_detail)

    for current_email_template in email_templates_list:
        # If the email template has a workflow name filter, skip if the submission doesn't have that workflow
        if not any(name in current_email_template.workflow_names for name in maybe_workflow_names):
            continue
        cli.logger.debug(f"Processing {current_email_template.config_name}")
        for submission_workflow in submission_detail["workflows"]:
            process_workflow(
                system_clients,
                gcloud_clients,
                submission_info,
                submission_workflow,
                current_email_template,
            )


def process_workflow(
        system_clients: SystemClients,
        gcloud_clients: GcloudClients,
        submission_info: SubmissionInfo,
        submission_workflow: dict[str, Any],
        current_email_template: EmailTemplate,
) -> None:
    workflow_id = submission_workflow["workflowId"]

    current_workflow_email = WorkflowEmail(current_email_template.config_name, workflow_id)
    if system_clients.data_store.has_workflow_email(current_workflow_email):
        cli.logger.debug(
            f"Already processed {current_email_template.config_name} in workflow {workflow_id}"
        )
        return

    cli.logger.debug(
        f"Processing {current_email_template.config_name} in workflow {workflow_id}"
    )
    workflow_metadata = gcloud_clients.terra.get_workflow(workflow_id)
    workflow_info = WorkflowInfo(
        submission_info.workspace_namespace,
        submission_info.workspace_name,
        submission_info.submission,
        submission_workflow,
        workflow_metadata,
    )
    process_workflow_email(
        system_clients,
        gcloud_clients,
        workflow_info,
        current_email_template,
        current_workflow_email,
    )


def process_workflow_email(
        system_clients: SystemClients,
        gcloud_clients: GcloudClients,
        workflow_info: WorkflowInfo,
        current_email_template: EmailTemplate,
        current_workflow_email: WorkflowEmail,
) -> None:
    output_paths = current_email_template.find_output_paths(workflow_info)
    missing_output_paths = [key for key, val in output_paths.items() if not val]
    if missing_output_paths:
        cli.logger.debug(
            f"Still missing outputs for {current_workflow_email.email_template}"
            f" in workflow {workflow_info.workflow_id}: {missing_output_paths}"
        )
        return
    output_data = {key: gcloud_clients.gcs.download_file(val) for key, val in output_paths.items()}
    email_message = current_email_template.make_message(workflow_info, output_data)
    send_email_message(system_clients, workflow_info, email_message, current_workflow_email)


def send_email_message(
        system_clients: SystemClients,
        workflow_info: WorkflowInfo,
        email_message: HtmlEmailMessage,
        current_workflow_email: WorkflowEmail,
) -> None:
    cli.logger.info(
        f"Sending {current_workflow_email.email_template} in workflow {workflow_info.workflow_id}"
    )
    if not email_message.email_from:
        email_message.email_from = system_clients.email_client.default_from_address()
    if not email_message.email_to:
        email_message.email_to = [workflow_info.submitter]
    system_clients.email_client.send_email(email_message)
    system_clients.data_store.add_workflow_email(current_workflow_email)


def run_with_retries(func: Callable[[], None], system_clients: SystemClients, recheck_seconds: Optional[int]) -> None:
    existing_error = False
    try:
        while True:
            try:
                func()
                if existing_error:
                    existing_error = False
                    cli.logger.info("Resuming normal operation")
            except KeyboardInterrupt:
                raise
            except Exception as exception:
                if not existing_error:
                    existing_error = True
                    cli.logger.exception(f"Error processing emails", exc_info=exception)
                    send_error_message(exception, system_clients)

            if recheck_seconds < 0:
                break
            cli.logger.debug(f"Sleeping for {recheck_seconds} seconds")
            time.sleep(recheck_seconds)
    except KeyboardInterrupt:
        pass


def send_error_message(exception: Exception, system_clients: SystemClients) -> None:
    if system_clients.email_client.errors_to_address():
        error_message = HtmlEmailMessage(
            email_from=system_clients.email_client.default_from_address(),
            email_to=[system_clients.email_client.errors_to_address()],
            subject="[INTERNAL_ERROR] : Trouble emailing outputs from Terra workflows",
            body=format_error_message(exception),
        )
        try:
            system_clients.email_client.send_email(error_message)
        except Exception as additional_exception:
            cli.logger.exception(f"Error sending error email", exc_info=additional_exception)


def format_error_message(exception: Exception) -> str:
    hr = ('<hr style="'
          ' border: 0;'
          ' height: 0;'
          ' border-top: 1px solid rgba(0, 0, 0, 0.1);'
          ' border-bottom: 1px solid rgba(255, 255, 255, 0.3);'
          '"/>')
    formatter = logging.Formatter()
    exception_string = formatter.formatException((type(exception), exception, exception.__traceback__))
    exception_html = html.escape(exception_string).replace("\n", "<br/>\n")
    return f"""
    <div>
        <div><table>
            <tr>
                <td style="text-align: right; font-weight: bold;">Status</td>
                <td><span style="color: #DF0101; font-weight: bold;">INTERNAL_ERROR</span></td>
            </tr>
        </table></div>
        {hr}
        <div>
            <div>{exception_html}</div>
        </div>
    </div>
    """


def main(args=None):
    parser = argparse.ArgumentParser(prog="email_outputs", description=__doc__)
    parser.add_argument("--log-level", "-l", default="INFO", choices=cli.dctLogLevel.keys(),
                        help="Set the logging level.  (default: %(default)s)")
    add_arguments(parser)

    if args is None:
        args = sys.argv[1:]

    options = parser.parse_args(args)
    cli.logger.setLevel(cli.dctLogLevel[options.log_level])
    return run(options)


if __name__ == "__main__":
    sys.exit(main())
