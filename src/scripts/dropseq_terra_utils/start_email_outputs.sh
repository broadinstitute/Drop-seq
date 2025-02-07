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

# Start the generic app and the loop monitor.
# The app specific code to run the app is contained in run_<app>.sh.
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
console_log="${root_dir}/log/console.log"

mkdir -p "${root_dir}/log"

force_start=0

function usage () {
    echo "USAGE: $prog_name [-f]" >&2
    echo  >&2
    echo "-f:    Force $app_name to start even if $app_name_lower.pid exists." >&2
}

while getopts ":fh" options; do
  case $options in
    f ) force_start=1;;
    h ) usage
          exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $((OPTIND - 1))

if [ -s "$app_pidfile" ]; then
  if [ "$force_start" == 1 ]; then
    echo "Removing existing $app_name_lower.pid ($(cat "$app_pidfile")) before launching." >&2
    rm "$app_pidfile"
  else
    echo "ERROR: $app_pidfile already exists. Is this instance already running? Stop the server or delete $app_pidfile before starting $app_name." >&2
    exit 1
  fi
fi

echo "Launching $app_name..."

set -m
nohup "${root_dir}/bin/loop_$app_name_lower.sh" >> "$console_log" 2>&1 &
set +m

loop_pid=$!
echo $loop_pid > "$loop_pidfile"
# Wait a few of seconds for launch
sleep 5
if [ ! -f "$app_pidfile" ]; then
	echo "ERROR: $app_name pid file not found: $app_pidfile." >&2
	exit 1
fi
app_pid=$(cat "$app_pidfile")
echo "$app_name launched with PID $app_pid and loop_$app_name_lower.sh PID $loop_pid..."

# Wait a little while for it to maybe fail before checking for process existence.
sleep 10

if ps -p "$app_pid" > /dev/null; then
  echo "$app_name appears to be running."
else
  echo "ERROR: $app_name process not found." >&2
  exit 1
fi
