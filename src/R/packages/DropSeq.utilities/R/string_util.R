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

#' Test if string ends with given suffix
#'
#' @param theString string to be tested
#' @param theExt suffix to be tested for
#' @return TRUE if theString ends with theExt
#' @export
strEndsWith<-function(theString,theExt) {
  return(substring(theString,1+nchar(theString)-nchar(theExt))==theExt)
}

#' Test if string starts with given prefix
#'
#' @param theString string to be tested
#' @param thePrefix prefix to be tested for
#' @return TRUE if theString starts with thePrefix
#' @export
strStartsWith = function(theString, thePrefix) {
  return(substr(theString, 1, nchar(thePrefix)) == thePrefix)
}

#' Remove a suffix from a string
#'
#' It is an error if the string does not end with the suffix.
#'
#' @param str the string to be truncated
#' @param suffix the string to be checked for and removed
#' @return the string with suffix removed
#' @export
removeSuffix = function(str, suffix) {
  if (!strEndsWith(str, suffix)) {
    stop(paste0("'", str, "' does not end with '", suffix, "'"))
  }
  return(substr(str, 1, nchar(str) - nchar(suffix)))
}

#' Remove a suffix from a string, if it is present.
#'
#' If the string does not end with the suffix, it is not truncated.
#'
#' @param str the string to be truncated
#' @param suffix the string to be checked for and removed
#' @return the string with suffix removed, if present, otherwise the original string
#' @export
maybeRemoveSuffix = function(str, suffix) {
  if (!strEndsWith(str, suffix)) {
    return(str)
  } else {
    return(substr(str, 1, nchar(str) - nchar(suffix)))
  }
}

#' Remove a suffix from a file and append a new one.
#' @param str String on which to operate.
#' @param oldSuffix suffix to be replaced.
#' @param newSuffix the replacement.
#' @return str with oldSuffix replaced by newSuffix
#' @export
replaceSuffix = function(str, oldSuffix, newSuffix) {
  return(paste0(removeSuffix(str, oldSuffix), newSuffix))
}

#' Split string into pieces of the given length
#' @param s String to be split
#' @param len desired length
#' @return vector of pieces.  Note that last piece may not be full length.
#' @export
splitStringByLength<-function(s, len) {
  if (nchar(s) <= len) return(s)
  starts <- seq(1, nchar(s), len)
  ends <- seq(len, nchar(s), len)
  if (length(ends) < length(starts)) {
    ends <- c(ends, nchar(s))
  }
  return(substring(s, starts, ends))
}

