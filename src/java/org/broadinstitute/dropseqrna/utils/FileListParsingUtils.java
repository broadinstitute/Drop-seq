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
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import picard.nio.PicardHtsPath;


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
        return expandPathList(IOUtil.filesToPaths(fileList)).stream().map(Path::toFile).collect(Collectors.toList());
    }

    /* Method returns the concatenated list of all the paths defined by the input paths
      Each input path is one of the following types:
      1. A path containing wildcard symbol(s) ("*" or "?"). E.g., /abc/file*.txt
        All the paths matching the wildcard pattern are added to the return file list.
      2. A path with the extension of '.bam_list' or '.file_list'.
        Each line of this path is treated as a path, and the corresponding path object is added to the return path list.
      3. Any other path is simply added to the returned path list.
     */
    public static List<PicardHtsPath> expandPicardHtsPathList(final List<PicardHtsPath> pathList) {
        return expandPathList(PicardHtsPath.toPaths(pathList))
                .stream()
                .map(PicardHtsPath::fromPath)
                .collect(Collectors.toList());
    }

    /* Method returns the concatenated list of all the paths defined by the input paths
      Each input path is one of the following types:
      1. A path containing wildcard symbol(s) ("*" or "?"). E.g., /abc/file*.txt
        All the paths matching the wildcard pattern are added to the return file list.
      2. A path with the extension of '.bam_list' or '.file_list'.
        Each line of this path is treated as a path, and the corresponding path object is added to the return path list.
      3. Any other path is simply added to the returned path list.
     */
    public static List<Path> expandPathList(final List<Path> pathList) {
        final List<Path> expandedPathList = new ArrayList<>();

        pathList.forEach(wildcardPath -> expandWildcardPath(wildcardPath).forEach(file -> {
            final String fileName = file.getFileName().toString();
            if (fileName.endsWith(".bam_list") || fileName.endsWith(".file_list"))
                expandedPathList.addAll(readPathList(file));
            else {
                IOUtil.assertFileIsReadable(file);
                expandedPathList.add(file);
            }
        }));

        return expandedPathList;
    }

    public static Collection<Path> expandWildcardPath(final Path wildcardPath) {
        try {
            Collection<Path> result;
            final String fileName = wildcardPath.getFileName().toString();
            if (fileName.contains("*") || fileName.contains("?")) {
                final WildcardFileFilter wildcardFileFilter =
                        WildcardFileFilter.builder().setWildcards(wildcardPath.getFileName().toString()).get();
                final Path parentPath =
                        wildcardPath.getParent() == null ? Paths.get(".") : wildcardPath.getParent();
                try (final DirectoryStream<Path> paths =
                             Files.newDirectoryStream(parentPath, wildcardFileFilter::matches)) {
                    result = Streams.stream(paths).collect(Collectors.toList());
                }
            } else {
                result = Collections.singleton(wildcardPath);
            }
            result.forEach(IOUtil::assertFileIsReadable);
            return result;
        } catch (IOException e) {
            throw new RuntimeIOException("Exception expanding wildcard path " + wildcardPath, e);
        }
    }

    /**
     * Implements the semantic in which a relative path in a file list is resolved relative to the canonical directory
     * containing the file list itself.
     */
    public static List<File> readFileList(final File fileList) {
        return readPathList(fileList.toPath()).stream().map(Path::toFile).collect(Collectors.toList());
    }

    /**
     * Implements the semantic in which a relative path in a file list is resolved relative to the canonical directory
     * containing the file list itself.
     */
    public static List<Path> readPathList(final Path pathList) {
        IOUtil.assertFileIsReadable(pathList);
        try {
            final Path canonicalDirectory = pathList.toRealPath().getParent();
            return Files.readAllLines(pathList)
                    .stream()
                    .map(canonicalDirectory::resolve)
                    .collect(Collectors.toList());
        } catch (IOException e) {
            throw new RuntimeIOException("Exception reading " + pathList, e);
        }
    }

    public static File resolveFilePath(final File fileListDirectory, final File filePath) {
        return fileListDirectory.toPath().resolve(filePath.toPath()).toFile();
    }
}
