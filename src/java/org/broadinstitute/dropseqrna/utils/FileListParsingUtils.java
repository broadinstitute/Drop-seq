/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;


public class FileListParsingUtils {

    /* Method returns the concatenated list of all the files defined by the input files
      Each input file is one of the following types:
      1. A file path containing wildcard symbol(s) ("*" or "?"). E.g., /abc/file*.txt
        All the files matching the wildcard pattern are added to the return file list.
      2. A file with the extension of '.bam_list' or '.file_list'.
        Each line of this file is treated as a file path, and the corresponding file object is added to the return file list.
      3. Any other file is simply added to the return file list.
     */
    public static List<File> expandFileList(final List<File> fileList) {
        final List<File> expandedFileList = new ArrayList<>();

        fileList.stream().forEach(wildcardFile -> expandWildcardFile(wildcardFile).stream().forEach(file -> {
            final String fileName = file.getName();
            if (fileName.endsWith(".bam_list") || fileName.endsWith(".file_list"))
                expandedFileList.addAll(readFileList(file));
            else {
            	IOUtil.assertFileIsReadable(file);
                expandedFileList.add(file);
            }
        }));

        return expandedFileList;
    }

    @SuppressWarnings("unchecked")
	public static Collection<File> expandWildcardFile(final File wildcardFile) {
    	Collection <File> result = new ArrayList<File>(); 
        final String fileName = wildcardFile.getName();
        if (fileName.contains("*") || fileName.contains("?")) {
            File parentFile = wildcardFile.getParentFile();
            if (parentFile == null) parentFile = new File(".");
            result = FileUtils.listFiles(parentFile, new WildcardFileFilter(wildcardFile.getName()), null);
        } else {
            result = Collections.singleton(wildcardFile);
        }
        result.stream().forEach(x -> IOUtil.assertFileIsReadable(x));
        return (result);
    }

    /**
     * Implements the semantic in which a relative path in a file list is resolved relative to the canonical directory
     * containing the file list itself.
     */
    public static List<File> readFileList(final File fileList) {
    	IOUtil.assertFileIsReadable(fileList);
        try {
            final File canonicalDirectory = fileList.getCanonicalFile().getParentFile();
            return Streams.stream((Iterable<String>) IOUtil.readLines(fileList)).
                    map(s -> resolveFilePath(canonicalDirectory, new File(s))).collect(Collectors.toList());
        } catch (IOException e) {
            throw new RuntimeIOException("Exception reading " + fileList, e);
        }
    }

    public static File resolveFilePath(final File fileListDirectory, final File filePath) {
        return fileListDirectory.toPath().resolve(filePath.toPath()).toFile();
    }
}
