#!/usr/bin/env bash
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

# Builds and installs the Drop-seq Python programs

# Uses uv to install the python environments.

# The multiple python environments should only depend on folders under $BASEDIR
# and not on any system python installations.

# Running pipx not only requires a system python, but will link created venvs to
# that system python unless we install yet another local python installation.
# We could use pyenv to create local pythons, but that installation depends on
# curl + git. uv only depends on curl so we're using that to install our venvs.

# If one can figure out the right incantation such that no files under $BASEDIR
# link to files outside of $BASEDIR, then we could use pipx install.

set -euo pipefail

uv_dir=/tmp/uv
envs_dir="$SRCDIR"/src/python
pythons_dir="$BASEDIR"/pythons
venv_dir=$BASEDIR/venv

export DEBIAN_FRONTEND=noninteractive
export UV_PYTHON_INSTALL_DIR=$pythons_dir
export UV_TOOL_DIR=$venv_dir
export UV_TOOL_BIN_DIR=$BASEDIR
export PATH=$PATH:$uv_dir:$BASEDIR

# Install prerequisites
echo "Installing prerequisites..."
apt-get -qq update && apt-get -qq install curl

# Install uv
curl -LsSf https://astral.sh/uv/install.sh | env UV_UNMANAGED_INSTALL=$uv_dir sh

# Install python scripts
for pyproject_toml in "$envs_dir"/*/pyproject.toml; do
  env_name=$(basename "$(dirname "$pyproject_toml")")
  echo "Building $env_name"
  uv tool install -q "$envs_dir/$env_name"
done
