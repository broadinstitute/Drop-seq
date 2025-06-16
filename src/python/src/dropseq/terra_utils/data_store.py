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

import os
from pathlib import Path
from typing import Optional, TypeVar, Generic, Callable

import pandas as pd
from pandas import Series

from dropseq.terra_utils.models import WorkflowEmail


def default_data_store_dir() -> Path:
    """
    Get the default data store directory.

    The data store directory is used to store data that persists between runs of the program.
    """
    data_home_dir = Path(os.getenv('XDG_STATE_HOME', Path.home() / '.local' / 'state'))
    data_store_dir = data_home_dir / 'dropseq_terra_utils'
    return data_store_dir


T = TypeVar('T')


class DataFile(Generic[T]):
    """
    A file that stores data that persists between runs of the program.
    The data is stored in a tab-separated file and loaded into memory when the file is accessed.
    """

    def __init__(
            self,
            file_path: Path,
            header: list[str],
            parse_row: Callable[[Series], T],
            to_row: Callable[[T], list[str]],
    ) -> None:
        self._file_path: Path = file_path
        self._header: list[str] = header
        self._parse_row: Callable[[Series], T] = parse_row
        self._to_row: Callable[[T], list[str]] = to_row
        self._cached_data: Optional[list[T]] = None

    def has_row(self, data: T) -> bool:
        """
        Check if the data is in the data file.
        """
        self.__ensure_data()
        return data in self._cached_data

    def add_row(self, data: T) -> None:
        """
        Add data to the data file.
        """
        self.__ensure_data()
        with open(self._file_path, "a") as file:
            file.write("\t".join(self._to_row(data)) + "\n")
        self._cached_data.append(data)

    def __ensure_data(self) -> None:
        """
        Ensure the existing data is loaded into memory.
        """
        if not self._file_path.exists():
            if not self._file_path.parent.exists():
                self._file_path.parent.mkdir(parents=True)
            with open(self._file_path, "w") as file:
                file.write("\t".join(self._header) + "\n")
        self._cached_data = (
            pd
            .read_csv(self._file_path, sep="\t")
            .apply(self._parse_row, axis=1)
            .tolist()
        )


class DataStore:
    """
    Used for storing data that persists between runs of the program.
    """

    def __init__(self, data_store_dir: Path):
        self.done_sumbissions_file: DataFile[str] = DataFile(
            data_store_dir / "done_submissions.tsv",
            ["submission_id"],
            lambda row: row["submission_id"],
            lambda data: [data],
        )
        self.workflow_emails_file: DataFile[WorkflowEmail] = DataFile(
            data_store_dir / "workflow_emails.tsv",
            ["email_template", "workflow_id"],
            lambda row: WorkflowEmail(row["email_template"], row["workflow_id"]),
            lambda data: [data.email_template, data.workflow_id],
        )

    def has_done_submission(self, submission_id: str) -> bool:
        return self.done_sumbissions_file.has_row(submission_id)

    def has_workflow_email(self, workflow_email: WorkflowEmail) -> bool:
        return self.workflow_emails_file.has_row(workflow_email)

    def add_done_submission(self, submission_id: str) -> None:
        self.done_sumbissions_file.add_row(submission_id)

    def add_workflow_email(self, workflow_email: WorkflowEmail) -> None:
        self.workflow_emails_file.add_row(workflow_email)
