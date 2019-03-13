#!/bin/bash
# MIT License
#
# Copyright 2017 Broad Institute
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

progname=__PROGNAME__
thisdir=`dirname ${BASH_SOURCE[0]}`
if [ -f $thisdir/loadDotKits.sh ]
then source $thisdir/loadDotKits.sh
fi

jar_deploy_dir=$thisdir/jar
verbose=0

USAGE=$(cat << EOF
USAGE: $0 [-m <jvm_heap_size>] [-v] program args...
       -m <jvm_heap_size> (default $xmx)
       -v echo final command line before executing

Program options:
EOF
)

function usage () {
    echo "$USAGE" >&2
    java -Xmx${xmx} -jar $jar_deploy_dir/dropseq.jar $progname -h
}

set -e

while getopts ":m:v" options; do
  case $options in
    m ) xmx=$OPTARG;;
    v ) verbose=1;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;

  esac
done
shift $(($OPTIND - 1))

broad_tmpdir_root=/broad/hptmp
if [ -z "$TMPDIR" -a -d $broad_tmpdir_root ]
then export TMPDIR=$broad_tmpdir_root/$USER
fi

if (( $verbose ))
then set -x
fi

java -Xmx${xmx} -XX:+UseParallelOldGC -XX:ParallelGCThreads=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+HeapDumpOnOutOfMemoryError -Djava.io.tmpdir=$TMPDIR -jar $jar_deploy_dir/dropseq.jar $progname $*