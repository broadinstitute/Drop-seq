# MIT License
#
# Copyright 2018 Broad Institute
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

# Shared definitions to be sourced by other scripts

thisdir=$(dirname $0)


function error_exit() {
    echo "ERROR: $@
    " >&2
    exit 1
}

if [[ -z "$verbose" ]]
then verbose=0
fi

if [[ -z "$ECHO" ]]
then ECHO=
fi

function check_invoke() {
    if (( $verbose ))
    then echo $@
    fi
    if $ECHO $@
    then :
    else error_exit "non-zero exit status " $? " executing $@"
    fi
}

picard_jar=$(find $thisdir -name picard\*.jar)

num_picard_jars=$(wc -w <<< "$picard_jar")

if (($num_picard_jars != 1))
then error_exit "Could not find one and only one picard.jar in deployment."
fi

function invoke_picard() {
    check_invoke java -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $picard_jar "$@"
}

function invoke_dropseq() {
    dropseq_program=$1
    shift
    check_invoke $thisdir/$dropseq_program "$@"
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

function check_TMPDIR() {
  if [[ -z "$TMPDIR" ]]
  then error_exit "TMPDIR environment variable must be set."
  fi
}