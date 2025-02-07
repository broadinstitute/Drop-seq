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

# Launch the generic app via run_<app>.sh and if it exits via signal (e.g. seg fault) restart it.
# NOTE: Start and stop the app using start_<app>.sh and stop_<app>.sh.
set -euo pipefail

cd "$(dirname "$0")/.."
prog_name=$(basename "$0")
root_dir=$(pwd)
app_name_lower=${prog_name}
app_name_lower=${app_name_lower%.sh}
app_name_lower=${app_name_lower#*_}
app_pidfile="${root_dir}/${app_name_lower}.pid"

exit_status=129

while (( exit_status > 128 )); do
  # shellcheck disable=SC1090
  source "${root_dir}/bin/run_${app_name_lower}.sh"

  app_pid=$!
  echo $app_pid > "$app_pidfile"

  # This waits until the process terminates, and set $? to the exit status of that program
  set +e
  wait $app_pid
  exit_status=$?
  set -e
done
