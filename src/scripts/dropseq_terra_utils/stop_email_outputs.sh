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

# Stop the generic app and the loop monitor.
set -euo pipefail

cd "$(dirname "$0")/.."
prog_name=$(basename "$0")
root_dir=$(pwd)
app_name_lower=${prog_name}
app_name_lower=${app_name_lower%.sh}
app_name_lower=${app_name_lower#*_}
app_name="$(tr '[:lower:]' '[:upper:]' <<< "${app_name_lower:0:1}")${app_name_lower:1}"
app_pidfile="${root_dir}/${app_name_lower}.pid"
loop_pidfile="${root_dir}/${app_name_lower}_loop.pid"

function process_exists {
    ps -p "$1" > /dev/null
}

if [ ! -s "$loop_pidfile" ]; then
  echo "WARNING: $loop_pidfile does not exist. Will not attempt to stop $app_name loop." >&2
else
  loop_pid=$(cat "$loop_pidfile")

  if process_exists "$loop_pid"; then
    kill -KILL "$loop_pid"
    echo "Sent stop signal to $app_name loop process $loop_pid."
  else
    echo "WARNING: $app_name loop process $loop_pid does not appear to exist. Continuing..." >&2
  fi
  rm "$loop_pidfile"
fi

if [ ! -s "$app_pidfile" ]; then
  echo "WARNING: $app_pidfile does not exist; assuming $app_name is already not running." >&2
  exit 0
fi
app_pid=$(cat "$app_pidfile")

if process_exists "$app_pid"; then
  kill -INT "$app_pid"
  echo "Sent stop signal to $app_name process $app_pid..."
else
  echo "$app_name process $app_pid does not appear to exist." >&2
  exit 1
fi

for _ in $(seq 10); do
  sleep 10
  process_exists "$app_pid" || {
    echo "$app_name process $app_pid has terminated."
    rm "$app_pidfile"
    exit 0
  }
  echo "Waiting for process to stop..."
done

echo "ERROR: $app_name process $app_pid did not terminate." >&2
exit 1
