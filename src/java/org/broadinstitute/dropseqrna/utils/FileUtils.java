/*
 * MIT License
 *
 * Copyright 2024 Broad Institute
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

import java.io.File;

public class FileUtils {
    /**
     * @return true if the given file ends with "." + periodlessExtension
     */
    public static boolean hasExtension(final File file,final String periodlessExtension) {
        return file.getName().endsWith("." + periodlessExtension);
    }

    /**
     *
     * @param file must end with "." + periodlessExtension
     * @param periodlessExtension
     * @return file object with the extension removed
     */
    public static File removeExtension(final File file,final String periodlessExtension) {
        if (!hasExtension(file, periodlessExtension)) {
            throw new IllegalArgumentException(String.format("%s does not have extension .%s", file, periodlessExtension));
        }
        final String filename = file.getName();
        return new File(file.getParentFile(), filename.substring(0, filename.length() - periodlessExtension.length() - 1));
    }

    /**
     * @return file with "." + periodlessExtension added
     */
    public static File addExension(final File file,final String periodlessExtension) {
        return new File(file.getParentFile(), file.getName() + "." + periodlessExtension);
    }

    public static File replaceExtension(final File file,final String periodlessExtensionToRemove, final String periodlessExensionToAdd) {
        return addExension(removeExtension(file, periodlessExtensionToRemove), periodlessExensionToAdd);
    }
}
