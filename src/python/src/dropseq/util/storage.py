#!/usr/bin/env python3
# MIT License
#
# Copyright 2026 Broad Institute
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
Read-only access to Google Cloud Storage (gs://) and the local filesystem through a common interface, so that
callers can work with either without localizing files.

Each store provides list_dir(url), read_text(url) and open_text(url).  Use default_store(url) to get the store
appropriate for a URL.
"""

import gzip
import io
import os

GCS_SCHEME = "gs://"


class GcsStore:
    """
    Read-only access to Google Cloud Storage.  The client is created on first use.
    """

    def __init__(self, client=None):
        self._client = client

    @property
    def client(self):
        if self._client is None:
            import google.cloud.storage
            self._client = google.cloud.storage.Client()
        return self._client

    @staticmethod
    def _split(url):
        bucket, _, name = url[len(GCS_SCHEME):].partition("/")
        return bucket, name

    def list_dir(self, url):
        """
        :return: (set of file basenames, sorted list of subdirectory names) directly under url.
        """
        bucket, prefix = self._split(url.rstrip("/") + "/")
        blobs = self.client.list_blobs(bucket, prefix=prefix, delimiter="/")
        files = {blob.name[len(prefix):] for blob in blobs if blob.name != prefix}
        # prefixes is only populated once the iterator has been consumed.
        subdirs = sorted(p[len(prefix):].rstrip("/") for p in blobs.prefixes)
        return files, subdirs

    def read_text(self, url):
        import google.cloud.storage
        return google.cloud.storage.Blob.from_string(url, self.client).download_as_text()

    def open_text(self, url):
        """
        :return: text file object for url, decompressed if url ends with .gz.
        """
        import google.cloud.storage
        return _text_stream(google.cloud.storage.Blob.from_string(url, self.client).open("rb"), url)


class LocalStore:
    """
    Access to the local filesystem, with the same interface as GcsStore.
    """

    def list_dir(self, url):
        if not os.path.isdir(url):
            return set(), []
        entries = os.listdir(url)
        files = {e for e in entries if os.path.isfile(os.path.join(url, e))}
        subdirs = sorted(e for e in entries if os.path.isdir(os.path.join(url, e)))
        return files, subdirs

    def read_text(self, url):
        with open(url) as f:
            return f.read()

    def open_text(self, url):
        """
        :return: text file object for url, decompressed if url ends with .gz.
        """
        return _text_stream(open(url, "rb"), url)


class _ClosingGzipFile(gzip.GzipFile):
    """
    GzipFile that also closes the file object it reads from.
    """

    def close(self):
        source = self.fileobj
        try:
            super().close()
        finally:
            if source is not None:
                source.close()


def _text_stream(binary, url):
    if url.endswith(".gz"):
        binary = _ClosingGzipFile(fileobj=binary)
    return io.TextIOWrapper(binary)


_stores = {}


def default_store(url):
    """
    :return: the store for url's storage type (GcsStore for gs:// URLs, otherwise LocalStore).  One store is
    shared per storage type, so a single GCS client is reused across calls.
    """
    is_gcs = url.startswith(GCS_SCHEME)
    if is_gcs not in _stores:
        _stores[is_gcs] = GcsStore() if is_gcs else LocalStore()
    return _stores[is_gcs]
