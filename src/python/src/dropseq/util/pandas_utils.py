#!/usr/bin/env python3
# MIT License
#
# Copyright 2022 Broad Institute
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
import sys

# Ensures that gzip files do not include the modified time in the header.
#
# See:
#  - https://docs.python.org/3/library/gzip.html#gzip.GzipFile.mtime
#  - https://github.com/openjdk/jdk8u/blob/jdk8u362-b00/jdk/src/share/classes/java/util/zip/GZIPOutputStream.java#L178-L194
gzip_compression = {'method': 'gzip', 'mtime': 0}


def get_compression(path):
    """Return the pandas compression type for a file."""
    if path and path.endswith('.gz'):
        return gzip_compression
    else:
        # Use the default compression arguments
        return 'infer'


def to_tsv(df, path_or_buf, **kwargs):
    """Write a dataframe to a tab separated file.

    Closes the path_or_buf when done except for sys.stdout or sys.stdout.buffer.
    """
    if path_or_buf == sys.stdout or (hasattr(sys.stdout, 'buffer') and path_or_buf == sys.stdout.buffer):
        # Write to stdout and leave the buffer open.
        df.to_csv(
            path_or_buf,
            sep='\t',
            **kwargs,
        )
    elif isinstance(path_or_buf, str):
        # Have pandas write to the path.
        df.to_csv(
            path_or_buf,
            sep='\t',
            compression=get_compression(path_or_buf),
            **kwargs,
        )
    else:
        # Close the existing buffer.
        path_or_buf.close()
        # Now tell pandas to write to the buffer's name.
        df.to_csv(
            path_or_buf.name,
            sep='\t',
            compression=get_compression(path_or_buf.name),
            **kwargs,
        )
