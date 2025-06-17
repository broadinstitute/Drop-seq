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
import shutil
import tempfile
import unittest
from pathlib import Path

from dropseq.terra_utils.models import WorkflowEmail
from dropseq.terra_utils.data_store import DataStore


class TestDataStore(unittest.TestCase):
    def setUp(self):
        self.tmpDir = Path(tempfile.mkdtemp(".tmp", "data_store."))

    def tearDown(self):
        shutil.rmtree(self.tmpDir)

    def test_workflow_email(self):
        existing_workflow_email = WorkflowEmail(
            email_template="existing_email_template",
            workflow_id="existing_workflow_id"
        )
        missing_workflow_email = WorkflowEmail(
            email_template="missing_email_template",
            workflow_id="missing_workflow_id"
        )

        ds = DataStore(self.tmpDir / "test_workflow_email")
        ds.add_workflow_email(existing_workflow_email)
        self.assertTrue(ds.has_workflow_email(existing_workflow_email))
        self.assertFalse(ds.has_workflow_email(missing_workflow_email))

        ds = DataStore(self.tmpDir / "test_workflow_email")
        self.assertTrue(ds.has_workflow_email(existing_workflow_email))
        self.assertFalse(ds.has_workflow_email(missing_workflow_email))

    def test_done_submission(self):
        ds = DataStore(self.tmpDir / "test_done_submission")
        ds.add_done_submission("existing_submission_id")
        self.assertTrue(ds.has_done_submission("existing_submission_id"))
        self.assertFalse(ds.has_done_submission("missing_submission_id"))

        ds = DataStore(self.tmpDir / "test_done_submission")
        self.assertTrue(ds.has_done_submission("existing_submission_id"))
        self.assertFalse(ds.has_done_submission("missing_submission_id"))


if __name__ == '__main__':
    unittest.main()
