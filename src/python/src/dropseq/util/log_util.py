#!/usr/bin/env python3
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

import logging
import time
import sys


def log_message(message, level=0, verbosity=0, file=sys.stderr, print_date=True, flush=False):
    if level <= verbosity:
        messages = [time.ctime(), message] if print_date else [message]
        print("\t".join(messages), file=file)
        if flush:
            file.flush()

# The below is a work in progress.  It needs to be encapsulated in a class or something to make it
# easier to use.
def configure_logger(name, level):
    # I cannot believe I need to do this to cause logger to write to stderr.
    logging.basicConfig(
        level=logging.INFO,               # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()] # StreamHandler writes to sys.stderr by default
    )
    return logging.getLogger(name)

def add_log_level_argument(parser, default="INFO"):
    """
    Add a log level argument to an argparse parser.

    Args:
        parser: The argparse parser to add the argument to.
        default: Default log level (default: "INFO").
    """
    parser.add_argument(
        "--log-level", "-l",
        choices=dctLogLevel.keys(),
        default=default,
        help="Set the logging level. (default: %(default)s)"
    )

def set_log_level_from_string(logger, logLevelStr):
    logger.setLevel(dctLogLevel[logLevelStr])

dctLogLevel = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL
}
