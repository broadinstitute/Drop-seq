# Miscellaneous python scripts for Drop-seq

## Installation

Requires python >= 3.12
```
pip install 'git+https://github.com/broadinstitute/Drop-seq.git#egg=dropseq&subdirectory=src/python'
```

## Usage

Wrapper scripts are generated for all the executables listed in `[project.scripts]` 
section of [pyproject.toml](pyproject.toml).  Each these will respond to `-h` or `--help`
with a usage message.