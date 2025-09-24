# Miscellaneous python scripts for Drop-seq

## Installation

Requires python >= 3.12
```
pip install 'git+https://github.com/broadinstitute/Drop-seq.git#egg=dropseq&subdirectory=src/python'
```
### Solving dependency issues when installing

This package is pure python, so should be installable on any platform with a working python 3.12+ 
installation. However, some of the dependencies may be difficult to install, typically because installation 
from sources  requires other tools on your system.  There are many python environment managers, so we 
won't provide manager-specific instructions, but here are some suggestions.

- For some dependencies, pre-built binary installations are available.  You can try to force pip to use
 those.  E.g., if you get errors installing numpy and pyarrow, you can add the argument 
`--only-binary numpy,pyarrow` to the pip install command.
- You can install the problematic dependencies separately, using the appropriate command for your
environment manager, and then install this package.
- You can disable dependency installation by adding `--no-deps` to the pip install command, and then
 install the problematic dependencies separately using the appropriate command for your environment
 manager.  You may also find installation instructions on the web sites for the dependencies. The list
 of packages dependencies can be found in the `dependencies` element of [pyproject.toml](pyproject.toml).
 Note also that depending on the scripts you want to run, you may not need all the dependencies.
- If an error message complains of the absence or wrong version of a tool like gcc or cmake, you can
  install the required version of that tool.

## Usage

Wrapper scripts are generated for all the executables listed in `[project.scripts]` 
section of [pyproject.toml](pyproject.toml).  Each these will respond to `-h` or `--help`
with a usage message.