#!/bin/sh
# MIT License
#
# Copyright 2023 Broad Institute
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


xmx=4g

progname=$(basename "$0")
thisdir=$(dirname "$0")
verbose=0

USAGE=$(cat << EOF
USAGE: $0 [-m <jvm_heap_size>] [-v] program args...

-m <jvm_heap_size> : Heap size to allocate for the JVM.  Default: $xmx.
-v                 : Echo final command line before executing.
-h                 : Print usage and exit.

Program options:
EOF
)

usage () {
    echo "$USAGE" >&2
    "$thisdir"/bin/dropseq "$progname" -h
}

set -e

while getopts ':m:vh' options; do
  case $options in
    m ) xmx=$OPTARG;;
    v ) verbose=1;;
    h ) usage
          exit 0;;
    \? ) break;; # exit option parsing at unrecognized option, so new-style Picard CLP args will survive.
    * ) usage
          exit 1;;

  esac
done
shift $((OPTIND - 1))

broad_tmpdir_root=/broad/hptmp
if [ -z "$TMPDIR" ] && [ -d $broad_tmpdir_root ]
then export TMPDIR="$broad_tmpdir_root/$USER"
fi


JAVA_OPTS="$JAVA_OPTS -Xmx$xmx"
if [ -n "$TMPDIR" ]
then JAVA_OPTS="$JAVA_OPTS -Djava.io.tmpdir=$TMPDIR"
fi

if [ "$verbose" -eq 1 ]
then set -x
fi

JAVA_OPTS="$JAVA_OPTS" "$thisdir"/bin/dropseq "$progname" "$@"
