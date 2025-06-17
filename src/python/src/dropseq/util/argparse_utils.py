#!/usr/bin/env python
# MIT License
# 
# Copyright 2019 Broad Institute
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
import sys
import argparse

def argparse_error(self, message):
    """Replacement method for ArgumentParser.error() the calls print_help rather than print_usage"""
    if message:
        self._print_message('ERROR: %s: %s\n\n' % (self.prog, message), sys.stderr)
    self.print_help(sys.stderr)
    if message:
        self._print_message('ERROR: %s: %s\n\n' % (self.prog, message), sys.stderr)
    self.exit(2)

class OptionFileArgumentParser(argparse.ArgumentParser):
    """
    Argument parser extension that handles picard-style options file.  NAME=VALUE options are parsed, one per line.
    NAME is lower-cased, underscores replaced with dashes, and preceded with dash or double dash depending on length
    of NAME.  These are then parsed by regular argument parser.
    """
    def __init__(self, *args, **kwargs):
        super(OptionFileArgumentParser, self).__init__(*args, fromfile_prefix_chars='@', **kwargs)
        # Capture the original file content
        self.argLines = []

    def convert_arg_line_to_args(self, arg_line):
        self.argLines.append(arg_line)
        fields = arg_line.split("=", 1)
        fields[0] = fields[0].lower().replace('_', '-')
        if len(fields[0]) == 1:
            fields[0] = '-' + fields[0]
        else:
            fields[0] = '--' + fields[0]
        if fields[0] in self._option_string_actions:
            return fields
        else:
            return []


