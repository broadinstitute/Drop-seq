#!/usr/bin/env python3
# MIT License
#
# Copyright 2024 Broad Institute
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
import anndata as ad
from typing import List, Set, Optional
import numpy as np
import pandas as pd
import logging
import re
logger = logging.getLogger(__name__)


def validate_expression_syntax(expr: str) -> None:
    """
    Check if a Python expression is syntactically valid.

    Parameters
    ----------
    expr : str
        The expression string to validate.

    Raises
    ------
    SyntaxError
        If the expression has invalid Python syntax.
    """
    try:
        compile(expr, "<filter_expr>", "eval")
    except SyntaxError as e:
        raise SyntaxError(f"Invalid expression syntax: {expr!r}\n{e}")
    
def check_filter_syntax(expr: str) -> bool:
    """
    Check if a filtering expression can be input directly to 
    evaluate_filters, i.e. if it is in the format:
    adata.obs['var'] <operator> 'value' 
    
    Parameters
    ----------
    expr : str
        The expression string to check.

    Returns
    -------
    bool
        True if the expression is in the correct format, False otherwise.
    """

    # Check if the expression contains 'adata.obs['
    if not re.search(r'adata\.obs\[[\'"]\w+[\'"]\]', expr):
        return False

    # Check for valid operators
    valid_operators = ['==', '!=', '<', '>', '<=', '>=']
    if not any(op in expr for op in valid_operators):
        return False

    # Check for valid values (strings or numbers)
    if not re.search(r'\'[^\']*\'|"[^"]*"', expr) and not re.search(r'\b\d+(\.\d+)?\b', expr):
        return False

    return True


def format_filter_syntax(expr: str) -> str:
    """
    Given a simple human-readable filtering expression, 
    return a string that can be fed to evaluate_filters directly, 
    in the format: adata.obs['var'] <operator> 'value'.

    Can handle simple human-readable expressions like:
    scpred_class == 'astrocyte'
    scpred_class == 'astrocyte' and scpred_class == 'oligodendrocyte'
    age >= 30 

    Parameters:
    ----------
    expr: str
        The human-readable filtering string to format.
        example: scpred_class == 'astrocyte'

    Returns:
    -------
    str
        The formatted filtering string to be used in evaluate_filters.
        example: adata.obs['scpred_class'] == 'astrocyte'
    """

    comp_pattern = r'(\b\w+\b)\s*([=!]=|[<>]=?)\s*(".*?"|\'.*?\'|\d+\.?\d*)'

    def replacer(match):
        var, op, val = match.groups()
        return f'(adata.obs["{var}"] {op} {val})'

    # replace all comparisons (e.g. ==, !=, <, >, <=, >=) with the correct syntax
    reformatted = re.sub(comp_pattern, replacer, expr)

    # replace logical operators with their Python equivalents
    reformatted = re.sub(r'\band\b', '&', reformatted)
    reformatted = re.sub(r'\bor\b', '|', reformatted)

    # ensure that reformatted expression is syntatically valid in Python
    validate_expression_syntax(reformatted)

    return reformatted


def evaluate_filters(adata: ad.AnnData, expressions: List[str]) -> ad.AnnData:

    """
    Apply filter expressions to an AnnData object and return a filtered copy.

    If no cells match, returns an empty AnnData object (same structure but zero observations).
    This avoids the need for None-checks downstream and ensures consistency in return type.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object (in-memory or backed).
    expressions : List[str]
        Python expressions evaluated in the context of the AnnData object.

    Returns
    -------
    AnnData
        A new AnnData object containing only the cells that match the filters.
        If no cells match, the returned object has the same structure but zero rows.

    Raises
    ------
    ValueError
        If any expression does not return a boolean array of the correct shape.
    """
    try:
        if not expressions:
            return adata

        keep = np.ones(adata.n_obs, dtype=bool)

        for expr in expressions:
            result = eval(expr, {}, {"adata": adata})

            if isinstance(result, pd.Series):
                result = result.to_numpy()

            if not isinstance(result, np.ndarray) or result.dtype != bool:
                raise ValueError(f"Expression did not return a boolean array: {expr!r}")
            if result.shape[0] != adata.n_obs:
                raise ValueError(f"Expression result length {result.shape[0]} does not match adata.n_obs")

            keep &= result

        if keep.all():
            return adata

        # If the AnnData object is backed, convert it to memory to avoid issues with slicing
        if adata.isbacked:
            adata = adata.to_memory()

        if not np.any(keep):
            return adata[[]].copy()

        return adata[keep]
    except Exception as e:
        raise ValueError(f"Failed to evaluate filters on AnnData object: {e}")

class _DummyAnnData:
    """Minimal stand-in so adata.obs resolves to obs_df in eval()."""
    def __init__(self, obs_df: pd.DataFrame):
        self.obs = obs_df

def filter_obs_dataframe(obs_df: pd.DataFrame, expressions: List[str]) -> pd.DataFrame:
    """
    Apply boolean filter expressions only to obs_df, supporting both
      - obs['col']
      - adata.obs['col']
    Raises ValueError on any bad expression, missing column, or non-boolean result.
    """
    if not expressions:
        return obs_df

    dummy = _DummyAnnData(obs_df)
    keep = np.ones(len(obs_df), dtype=bool)

    for expr in expressions:
        try:
            result = eval(expr, {}, {"obs": obs_df, "adata": dummy})
        except KeyError as e:
            # missing column in obs_df
            raise ValueError(f"Column {e.args[0]!r} not found in obs when evaluating filter {expr!r}")
        except Exception as e:
            # syntax errors, name errors, etc.
            raise ValueError(f"Error evaluating filter {expr!r}: {e}")

        if isinstance(result, pd.Series):
            result = result.to_numpy()

        if not (isinstance(result, np.ndarray)
                and result.dtype == bool
                and result.shape[0] == len(obs_df)):
            raise ValueError(
                f"Filter {expr!r} must yield a boolean array of length {len(obs_df)}"
            )

        keep &= result

    if not keep.any():
        return obs_df.iloc[0:0]
    return obs_df.loc[keep]

def check_and_reformat_filter(expr: str) -> str:
    """
    Check if the expression is in the correct format for evaluation
    (eg. adata.obs['var'] <operator> 'value').
    
    If not, reformat it to the correct format.

    Parameters
    ----------
    expr : str
        The expression string to check and reformat.
    
    Returns
    -------
    str
        The reformatted expression string.
    """

    if check_filter_syntax(expr):
        return expr
    else:
        return format_filter_syntax(expr)

def evaluate_any_filters(adata: ad.AnnData, expressions: List[str]) -> ad.AnnData:
    """
    Evaluate ANY filter expressions on an AnnData object. 
    
    These expressions can be in EITHER the human-readable format (e.g. scpred_class == 'astrocyte')
    or the format (e.g. adata.obs['scpred_class'] == 'astrocyte'). If human-readable, the 
    expression is reformatted to the latter before evaluation.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object (in-memory or backed).
    
    expressions : List[str]
        List of filter expressions to evaluate.
        
    Returns
    -------
    AnnData
        A new AnnData object containing only the cells that match the filters.
        If no cells match, the returned object has the same structure but zero rows.
    """

    # first check if the expressions are syntactically valid in Python
    for expr in expressions:
        validate_expression_syntax(expr)
    
    # reformat expressions that are not in the correct format
    expressions= [check_and_reformat_filter(expr) for expr in expressions]

    # apply filters to anndata object
    return evaluate_filters(adata, expressions)
