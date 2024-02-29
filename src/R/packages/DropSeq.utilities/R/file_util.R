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

#' Read a data.table using data.table fread()
#'
#' If inFile ends with .gz, it is gunzipped in a way that is robust to running out of disk space.
#' Environment variable DROPSEQ_TMPDIR should be set to a directory that has plenty of room and is writable.
#'
#' @param inFile path of file to be read
#' @param comment_regexp If defined, filter input through "egrep -v comment_regexp"
#' @param ... passed through to data.table::fread
#' @return a data.table
#' @import data.table
#' @export
fastRead<-function(inFile, comment_regexp=NA, ...) {
  if (!file.exists(inFile)) {
    stop(paste0(inFile, " does not exist"))
  }
  # fread writes process output to /dev/shm if it exists, which is OK so long as it is big enough.
  # If that fails, fall back to function that uses an explicit TMPDIR
  if (length(grep(".gz", inFile))==1) {
    unzip_command = paste("gunzip","-c",inFile,sep=" ")
    if (!is.na(comment_regexp)) {
      unzip_command = paste0(unzip_command, sprintf(" | egrep -v '%s'", comment_regexp))
    }
    # data.table 1.11.8 emits a warning if command string isn't passed via cmd=, but v 1.10.4 (current dotkit version)
    # doesn't recognized cmd=
    if (version_above("data.table", "1.11.8")) {
      return(tryCatch(data.table::fread(cmd=unzip_command,data.table=T, ...), error=function(e) fastReadBigGz(inFile, comment_regexp=comment_regexp, ...)))
    } else {
      return(tryCatch(data.table::fread(unzip_command,data.table=T, ...), error=function(e) fastReadBigGz(inFile, comment_regexp=comment_regexp, ...)))
    }

  } else {
    input = inFile
    if (!is.na(comment_regexp)) {
      input = sprintf("egrep -v '%s' %s", comment_regexp, inFile)
      a=data.table::fread(input, data.table=T, ...)
    } else {
      a=data.table::fread(input, data.table=T, ...)
    }
    return(a)
  }
}

# internal function
fastReadBigGz<-function(inFile, ...) {
  # Grid Engine overwrites TMPDIR, so support our own environment variable that won't get clobbered.
  tmpDir = Sys.getenv('DROPSEQ_TMPDIR')
  if (nchar(tmpDir) == 0) {
    tmpDir = Sys.getenv('TMPDIR')
  }
  if (nchar(tmpDir) == 0) {
    tmpDir = tempdir()
  }
  t=tempfile(tmpdir=tmpDir, fileext=".tsv")
  on.exit(unlink(t), add = TRUE)
  cmd=paste("gunzip -c", inFile, ">", t)
  retval = system(cmd,intern=F) != 0
  if (retval != 0) {
    stop(paste0(cmd, " failed with status ", retval))
  }
  return(fastRead(t, ...))
}

#' Read a tabular DGE, with gene rows and columns that typically are cell barcodes
#' @param file to be read, optionally gzipped
#' @param decreasing_order_by_size If true, columns are ordered in decreasing order by colsum
#' @return a data.frame with genes as rownames, no GENE column.
#' @export
#' @import data.table
read_dge_gz<-function(file, decreasing_order_by_size=TRUE) {
  dge = fastRead(file, comment_regexp = '^#')
  setkey(dge,GENE)
  gene_names = dge$GENE
  GENE <- NULL # Silence R CMD check warning https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  dge[, GENE:=NULL]
  rownames(dge) = gene_names

  if (decreasing_order_by_size) {
    # order cells by expression

    # I don't understand data.table, but this returns a data.table with one row rather than
    # a vector, so convert to list to avoid warnings
    colsums.full.data=unlist(dge[,lapply(.SD,sum,na.rm=TRUE),.SDcols=colnames(dge)])
    order.data=names(colsums.full.data)[order(colsums.full.data,decreasing=T)]
    full.data.cells=dge[,order.data,with=F]
  } else {
    full.data.cells = dge
  }
  full.data.cells = as.data.frame(full.data.cells)
  rownames(full.data.cells) = gene_names
  return(full.data.cells)
}

#' Creates a file connection with the given open mode.
#' @param file If ends with ".gz", a gzfile() is created; else a regular file() connection.
#' @param open mode in which file is opened.  Default: "rb"
#' @export
open_conn = function(file, open="") {
  if (strEndsWith(file, ".gz")) {
    # work around bug in gzfile:
    # https://stackoverflow.com/questions/45665496/how-would-one-readlines-from-a-gzip-file-in-r
    if (nchar(open) == 0) {
      open = "rb"
    }
    if (!strEndsWith(open, "b")) {
      open=paste0(open, "b")
    }
    return(gzcon(file(file, open=open)))
  } else {
    return(file(file, open=open))
  }
}

#' Write a tab-separated file with row and column labels, including a label for the rownames.
#' Additional args are passed to write.table
#' @param x the data frame to be written.
#' @param file where to write the data.  If suffix is ".gz", the file is gzipped.
#' @param rowname_label the column header for the rownames
#' @param quote passed to write.table
#' @param sep passed to write.table
#' @param col.names passed to write.table
#' @param ... passed to write.table
#' @export
write_table_with_label_for_rownames = function(x, file, rowname_label, quote=FALSE, sep="\t", col.names=TRUE, ...) {
  conn = open_conn(file=file, open="w")
  cat(paste0(rowname_label, '\t'), file=conn)
  utils::write.table(x, conn, col.names=col.names, row.names=TRUE, quote=quote, sep=sep, ...)
  close(conn)
}

#' Read a tab-separated file, drop the first column and use the values of the first column as rownames.
#' @param file If suffix is ".gz" it is assumed the file is gzipped.
#' @param ... passed to data.table::fread()
#' @export
#' @import data.table
read_table_with_label_for_rownames=function(file, ...) {
  ret = fastRead(file, ...)
  rowname_column = colnames(ret)[1]
  labels = ret[, get(rowname_column)]
  ret[, (rowname_column):=NULL]
  ret = as.data.frame(ret)
  rownames(ret) = labels
  return(ret)
}

#' Make sure the argument has a trailing slash.
#' @param directory path on which to operate
#' @return the input argument, with a slash appended if it doesn't already end with a slash
#' @export
ensure_trailing_slash=function(directory) {
  if (!strEndsWith(directory, "/")) {
    directory = paste0(directory, "/")
  }
  return(directory)
}

#' Check if a directory exists, and create it if not.  Errors if cannot be created or exists but is not a directory.
#' @param directory to be created
#' @param recursive if FALSE, the parent directory must exist.
#' @return the directory, always with a trailing slash.  Raises an error if cannot be created
#' @export
ensure_directory_exists = function(directory, recursive=TRUE) {
  directory = ensure_trailing_slash(directory)
  if (dir.exists(directory)) {
    return(directory);
  } else if (file.exists(directory)) {
    stop(paste0("Cannot create directory because there is an existing file with the same name: ", directory))
  } else if (!dir.create(directory, recursive=recursive)) {
    stop(paste0("Error creating directory: ", directory))
  } else {
    return(directory)
  }
}

#' Count the number of lines at the beginning of the file that start with @ or #
#' @param path text file to be read, possibly gzipped
#' @return the number of lines at the beginning of the file that start with '@' or '#'
#' @export
countHeaderLines = function(path) {
  conn = file(path, "r")
  on.exit(close(conn), add = TRUE)
  headerLength = 0
  while(strStartsWith(line<-(readLines(conn, n=1)[[1]]), '@') || strStartsWith(line, '#')) {
    headerLength = headerLength + 1
  }
  return(headerLength)
}

dge_header_attribute_name='Drop-seq.dge.header'

version_above <- function(pkg, than) {
  as.logical(utils::compareVersion(as.character(utils::packageVersion(pkg)), than))
}

#' Read a tabular DGE, that may have comment/header lines starting with '@' or '#'
#' @param path file to be read.  If has suffix ".gz" it is assumed to be gzipped.
#' @return a data.table with the tabular contents of the file, and header/comment lines attached
#'         as an attribute named 'Drop-seq.dge.header'
#' @import data.table
#' @export
read_dge_txt = function(path) {
  headerLength = countHeaderLines(path)
  dge = fastRead(path, skip=headerLength)
  if (headerLength > 0) {
    header = readLines(path, n=headerLength)
    attr(dge, dge_header_attribute_name) = header
  }
  return(dge)
}

#' Determine if a file is in matrix market file format.
#'
#' @param file The file to test
#'
#' @return true if the file is in matrix market format
#' @export
isMatrixMarket<-function(file) {
  conn=file(file, "r")
  line=readLines(con=conn, n=1, ok=FALSE)
  close(conn)
  return(strStartsWith(line, "%%MatrixMarket matrix coordinate"))
}

#' read the first fiew lines of a MatrixMarket file to get the header line.
#' @param file possibly gzipped
#' @return list(numRows, numCols, numPopulatedElements)
#' @export
getSparseMatrixDimensions=function(file) {
  conn=open_conn(file, "r")
  line = "%"
  while(strStartsWith(line, '%')) {
    line=readLines(con=conn, n=1, ok=FALSE)
  }
  fields = strsplit(line, "\t", fixed=TRUE)[[1]]
  return(list(numRows=as.numeric(fields[1]), numCols=as.numeric(fields[2]), numPopulatedElements=as.numeric(fields[3])))
}

#' Read enough of a file to get the number of rows and columns.
#' @param file possibly gzipped
#' @return list(numRows, numCols)
#' @export
getDenseDgeDimensions=function(file) {
  numCols = ncol(fastRead(file, nrows=1)) - 1
  # start at -1 because first line is a header
  numRows = -1
  conn = open_conn(file)
  while (TRUE) {
    l = readLines(conn, ok=TRUE, n=1)
    if (length(l) == 0) {
      break
    }
    if (!strStartsWith(l, "#")) {
      numRows = numRows + 1
    }
  }
  close(conn)
  return(list(numRows=numRows, numCols=numCols))
}

#' Determine size of sparse or dense DGE
#' @param file possibly gzipped
#' @return list(numRows, numCols, numPopulatedElements)
#' @export
getDgeDimensions<-function(file) {
  if (isMatrixMarket(file)) {
    return(getSparseMatrixDimensions(file))
  } else {
    ret = getDenseDgeDimensions(file)
    ret$numPopulatedElements = ret$numRows * as.numeric(ret$numCols)
    return(ret)
  }
}

#' Strings that appear in MatrixMarket files.
#' @export
MatrixMarketConstants = list(
  MM_COMMENT_LINE_START = "%",
  MM_STRUCTURED_COMMENT_LINE_START = "%%",
  MM_MATRIX_TYPE_GENERAL="general",
  MM_HEADER_LIST_SEPARATOR = "\t",
  real = "real",
  integer = "integer",
  real_element_format = "%d\t%d\t%.8g",
  integer_element_format = "%d\t%d\t%d"
)
MatrixMarketConstants$MM_HEADER_START = paste0(MatrixMarketConstants$MM_STRUCTURED_COMMENT_LINE_START, "MatrixMarket matrix coordinate")
