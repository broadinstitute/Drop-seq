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
progname=$(basename "$0")
main_class=org.broadinstitute.dropseqrna.cmdline.DropSeqMain
usage () {
    cat >&2 <<EOF
USAGE: $progname -t <template-file> -c <classpath> -d <output-directory> -m <main-class> [main-class-args...]
Create wrapper scripts for Java command-line programs

-t <template-file>  : File to be copied to make the wrapper.  Required.
-c <classpath>      : Classpath for running the CLP lister.  Required.
-d <output-directory>  : Where to write the wrappers.  Required.
-m <main-class>     : Class to invoke to list CLPs.  Default: $main_class.
-h                  : Print usage and exit.
[main-class-args]   : Passed to main-class
EOF
}

error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

check_set() {
    value=$1
    name=$2
    flag=$3

    if [ -z "$value" ]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e

while getopts ':t:c:d:m:h' options; do
  case $options in
    t ) template=$OPTARG;;
    c ) classpath=$OPTARG;;
    d ) outdir=$OPTARG;;
    m ) main_class=$OPTARG;;
    h ) usage
          exit 0;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $((OPTIND - 1))

check_set "$outdir" 'output directory' '-d'
check_set "$template" 'wrapper script template' '-t'
check_set "$classpath" 'classpath' '-c'

main_class_args="$*"
mkdir -p "$outdir"
clp_names=$(java -cp "$classpath" "$main_class" "$main_class_args" || echo > /dev/null)
if [ -z "$clp_names" ]
then echo 'There was a problem getting CLP names.'
     exit 1
fi
for clp_name in $clp_names
do cp -p "$template" "$outdir/$clp_name"
done

echo 'Created wrapper scripts'
