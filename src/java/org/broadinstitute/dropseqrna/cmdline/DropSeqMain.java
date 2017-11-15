/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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
package org.broadinstitute.dropseqrna.cmdline;

import picard.cmdline.PicardCommandLine;

import java.util.ArrayList;
import java.util.List;

public class DropSeqMain extends PicardCommandLine {
    /** The name of this unified command line program **/
    private final static String COMMAND_LINE_NAME = DropSeqMain.class.getSimpleName();

    /** The packages we wish to include in our command line **/
    protected static List<String> getPackageList() {
        final List<String> packageList = new ArrayList<String>();
        packageList.add("org.broadinstitute.dropseqrna");
        return packageList;
    }
    public static void main(final String[] args) {
        System.exit(new DropSeqMain().instanceMain(args, getPackageList(), COMMAND_LINE_NAME));
    }

}
