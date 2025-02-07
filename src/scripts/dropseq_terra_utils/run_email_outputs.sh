#!/bin/bash

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

# Run the specified application with custom commands.
# Sourced and monitored by the loop_<app>.sh sibling script.
# Start and stop the app using start_<app>.sh and stop_<app>.sh.

set -euo pipefail

cd "$(dirname "$0")/.."
root_dir=$(pwd)
dropseq_dir=/broad/mccarroll/software/dropseq/priv
project_metadata_yaml="${root_dir}/conf/project_metadata.yaml"

# Using intermediate bash scripts interferes with the SIGINT sent by stop_email_outputs.sh.
# Do all the setup here and then directly launch the python process.
set +u
source $dropseq_dir/configDropSeqRNAEnvironment.bash
set -u
source activate base
source activate $dropseq_dir/conda/envs/dropseq_terra_utils
export PYTHONPATH=$dropseq_dir/pythonuserbase/lib/python3.12/site-packages

# enable job control so that kill -INT can be sent to background process
set -m
set -x
python -m dropseq_terra_utils.email_outputs --project-metadata "$project_metadata_yaml" &
set +x
set +m
