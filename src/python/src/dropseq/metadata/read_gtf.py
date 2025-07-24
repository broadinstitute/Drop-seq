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

# Based on https://github.com/openvax/gtfparse

import pandas as pd

class GtfRequiredColNames:
    CHROMOSOME = 'CHROMOSOME'
    SOURCE = 'SOURCE'
    FEATURE= 'FEATURE'
    START = 'START'
    END = 'END'
    SCORE = 'SCORE'
    STRAND = 'STRAND'
    FRAME = 'FRAME'
    ATTRIBUTE = 'ATTRIBUTE'
    all_colnames = [
        CHROMOSOME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE
    ]

def fix_attribute_field(attribute):
    """
    Fix the attribute field in a GTF file by replacing quotes and semicolons.

    Args:
        attribute (str): The attribute field from a GTF file.

    Returns:
        str: The fixed attribute field.
    """
    # Replace quotes with empty string and semicolons with commas
    return attribute.replace('"', "'").replace(';\"', '\"').replace(";-", "-")

def parse_attribute_to_dict(attribute):
    """
    Parse the attribute field of a GTF file into a dictionary.

    Args:
        attribute (str): The attribute field from a GTF file.

    Returns:
        dict: A dictionary with parsed attributes.
    """
    # Fix the attribute field
    fixed_attribute = fix_attribute_field(attribute)

    # Split by semicolon and strip whitespace
    pairs = [pair.strip() for pair in fixed_attribute.split(';') if pair.strip()]

    # Create a dictionary from the pairs
    return {key_value.split(' ')[0]: key_value.split(' ')[1].strip("'") for key_value in pairs}

def read_gtf(gtf_path, filterFunc = lambda df: df):
    """
    Read a GTF file and return a DataFrame with gene information.

    Args:
        gtf_path (str): Path to the GTF file.
        filterFunc (function): A function to filter the DataFrame after reading. May be used to reduce the number of
        rows before parsing the attributes, or eliminate ATTRIBUTE column, to speed loading. Defaults to a no-op function.

    Returns:
        pd.DataFrame: DataFrame containing gene information.
    """

    gtf = pd.read_csv(gtf_path, sep='\t', header=0, index_col=None, names=GtfRequiredColNames.all_colnames, comment='#')
    gtf = filterFunc(gtf)

    if GtfRequiredColNames.ATTRIBUTE in gtf.columns:
        parsedAttributes = gtf[GtfRequiredColNames.ATTRIBUTE].apply(parse_attribute_to_dict)
        # Drop the original ATTRIBUTE column
        gtf.drop(columns=[GtfRequiredColNames.ATTRIBUTE], inplace=True)
        # Create a DataFrame from the parsed attributes
        attributes_df = pd.DataFrame(parsedAttributes.tolist())
        gtf = pd.concat([gtf.reset_index(drop=True), attributes_df.reset_index(drop=True)], axis=1)
    return gtf



