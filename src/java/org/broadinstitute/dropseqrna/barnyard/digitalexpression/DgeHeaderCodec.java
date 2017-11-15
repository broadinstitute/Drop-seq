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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.util.*;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;

public class DgeHeaderCodec {

    private static final Log log = Log.getInstance(DgeHeaderCodec.class);

    private static char RECORD_START = '@';
    private static String KV_SEPARATOR = ":";
    private static String FIELD_SEPARATOR = "\t";
    private static String DGE_RECORD_LABEL = "DGE";
    private static String LIBRARY_RECORD_LABEL = "LIBRARY";
    private static String COMMAND_RECORD_LABEL = "COMMAND";
    private enum DgeRecordTag {VERSION, EXPRESSION_FORMAT}
    private enum LibraryRecordTag {INPUT, INPUT_DGE, REFERENCE, UEI, PREFIX}
    private enum CommandRecordTag {CL};

    public void encode(final Writer writer, final DgeHeader header) {
        writeLine(writer, buildFirstLine(header));
        for (final String prefix : new IterableAdapter<>(header.iterateLibraries())) {
            writeLine(writer, buildLibraryLine(header.getLibrary(prefix)));
        }
        for (final DgeHeaderCommand command: new IterableAdapter<>(header.iterateCommands())) {
            writeLine(writer, buildCommandLine(command));
        }
    }

    public void encode(final File file, final DgeHeader header) {
        final Writer writer = IOUtil.openFileForBufferedWriting(file);
        encode(writer, header);
        try {
            writer.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + file.getAbsolutePath(), e);
        }
    }

    private void writeLine(final Writer writer, final String line) {
        try {
            writer.write(line);
            writer.write('\n');
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    private String buildFirstLine(final DgeHeader header) {
        final OutputRecordBuilder dgeRecord = new OutputRecordBuilder(DGE_RECORD_LABEL);
        dgeRecord.addFieldIfNotNull(DgeRecordTag.VERSION, header.getVersion());
        dgeRecord.addFieldIfNotNull(DgeRecordTag.EXPRESSION_FORMAT, header.getExpressionFormat());
        return dgeRecord.build();
    }

    private String buildLibraryLine(final DgeHeaderLibrary library) {
        final OutputRecordBuilder dgeRecord = new OutputRecordBuilder(LIBRARY_RECORD_LABEL);
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.INPUT, library.getInput());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.INPUT_DGE, library.getInputDge());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.REFERENCE, library.getReference());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.UEI, library.getUei());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.PREFIX, library.getPrefix());
        for (final Map.Entry<String, String> entry: new IterableAdapter<>(library.iterateOtherTags())) {
            dgeRecord.addFieldIfNotNull(entry.getKey(), entry.getValue());
        }
        return dgeRecord.build();
    }

    private String buildCommandLine(final DgeHeaderCommand command) {
        final OutputRecordBuilder dgeRecord = new OutputRecordBuilder(COMMAND_RECORD_LABEL);
        dgeRecord.addFieldIfNotNull(CommandRecordTag.CL, command.getCommandLine());
        return dgeRecord.build();
    }

    private static class OutputRecordBuilder {
        private final ArrayList<String> fields = new ArrayList<>();

        public OutputRecordBuilder(String recordType) {
            fields.add(RECORD_START + recordType);
        }

        public void addFieldIfNotNull(final Object key, final Object value) {
            if (value != null) {
                fields.add(key + KV_SEPARATOR + value);
            }
        }

        public String build() {
            return StringUtil.join(FIELD_SEPARATOR, fields);
        }
    }


    /**
     *
     * @param reader When method returns, reader is positioned at start of first line after header
     */
    public DgeHeader decode(final BufferedReader reader, final String inputName) {
        if (!reader.markSupported()) {
            throw new IllegalArgumentException("reader.markSupported == false");
        }
        return decode(new DgeBufferedReader(reader), inputName);
    }
    /**
     *
     * @param inputStream When method returns, input stream is positioned at start of first line after header
     */
    public DgeHeader decode(final BufferedInputStream inputStream, final String inputName) {
        if (!inputStream.markSupported()) {
            throw new IllegalArgumentException("reader.markSupported == false");
        }
        return decode(new DgeBufferedInputStream(inputStream), inputName);
    }



    private DgeHeader decode(final DgeHeaderLineReader reader, final String inputName) {
        try {
            final DgeHeader ret = new DgeHeader();
            ret.setExpressionFormat(DgeHeader.ExpressionFormat.unknown);
            ret.setVersion(null);
            boolean first = true;
            String line;
            while ((line = reader.readHeaderLine()) != null) {
                if (line.startsWith(DGE_RECORD_LABEL + FIELD_SEPARATOR)) {
                    if (first) {
                        parseFirstLine(ret, line, inputName);
                    } else {
                        log.warn("DGE header line seen after first line in " + inputName);
                    }
                } else if (line.startsWith(LIBRARY_RECORD_LABEL + FIELD_SEPARATOR)) {
                    DgeHeaderLibrary library = parseLibrary(line, inputName);
                    if (library != null) {
                        ret.addLibrary(library);
                    }
                } else if (line.startsWith(COMMAND_RECORD_LABEL + FIELD_SEPARATOR)) {
                    DgeHeaderCommand command = parseCommand(line, inputName);
                    if (command != null) {
                        ret.addCommand(command);
                    }
                } else {
                    log.warn("Strange DGE header line seen in " + inputName + ";" + line);
                }
                first = false;
            }
            return ret;
        } catch (IOException e) {
            throw new RuntimeException("Problem parsing " + inputName, e);
        }
    }

    private interface DgeHeaderLineReader {
        /**
         * @return next header line, without leading '@', or null if no more header lines
         */
        String readHeaderLine() throws IOException;
    }

    private class DgeBufferedReader implements DgeHeaderLineReader {
        private final BufferedReader reader;

        public DgeBufferedReader(BufferedReader reader) {
            this.reader = reader;
        }

        @Override
        public String readHeaderLine() throws IOException {
            reader.mark(1);
            int c = reader.read();
            if (c == RECORD_START) {
                return reader.readLine();
            } else {
                reader.reset();
                return null;
            }
        }
    }

    private class DgeBufferedInputStream implements DgeHeaderLineReader {
        private final BufferedInputStream inputStream;
        // Could be local to readHeaderLine, but time is saved by allocating once.
        private char[] lineBuffer = new char[10000];;

        public DgeBufferedInputStream(BufferedInputStream inputStream) {
            this.inputStream = inputStream;
        }
        @Override
        public String readHeaderLine() throws IOException {
            inputStream.mark(1);
            int c = inputStream.read();
            if (c == RECORD_START) {
                return readToEndOfLine();
            } else {
                inputStream.reset();
                return null;
            }
        }

        private static final int BUFFER_OVERFLOW_INCREASE_FACTOR = 2;
        private static final byte LINEFEED = (byte) ('\n' & 0xff);
        private static final byte CARRIAGE_RETURN = (byte) ('\r' & 0xff);

        // Cribbed from Picard AsciiLineReader
        private String readToEndOfLine() throws IOException {
            int linePosition = 0;

            while (true) {
                final int b = inputStream.read();

                if (b == -1) {
                    // eof reached.  Return the last line, or null if this is a new line
                    if (linePosition > 0) {
                        return new String(lineBuffer, 0, linePosition);
                    } else {
                        return null;
                    }
                }

                final char c = (char) (b & 0xFF);
                if (c == LINEFEED || c == CARRIAGE_RETURN) {
                    if (c == CARRIAGE_RETURN && peek() == LINEFEED) {
                        inputStream.read(); // <= skip the trailing \n in case of \r\n termination
                    }

                    return new String(lineBuffer, 0, linePosition);
                } else {
                    // Expand line buffer size if neccessary.  Reserve at least 2 characters
                    // for potential line-terminators in return string

                    if (linePosition > (lineBuffer.length - 3)) {
                        final char[] temp = new char[BUFFER_OVERFLOW_INCREASE_FACTOR * lineBuffer.length];
                        System.arraycopy(lineBuffer, 0, temp, 0, lineBuffer.length);
                        lineBuffer = temp;
                    }

                    lineBuffer[linePosition++] = c;
                }
            }

        }

        private int peek() throws IOException {
            inputStream.mark(1);
            try {
                return inputStream.read();
            } finally {
                inputStream.reset();
            }
        }
    }

    /**
     *
     * Checks the next character.  If '@', that char is consumed.  If not, the char read is pushed back into the reader.
     * @return true if next character is '@'
     */
    private boolean isHeaderLine(final BufferedReader reader) throws IOException {
        reader.mark(1);
        int c = reader.read();
        if (c == RECORD_START) {
            return true;
        } else {
            reader.reset();
            return false;
        }
    }

    private void parseFirstLine(final DgeHeader header, final String line, final String inputName) {
        final LinkedHashMap<String, String> fields = parseLine(line, DGE_RECORD_LABEL, inputName);
        for (final Map.Entry<String, String> entry: fields.entrySet()) {
            final String key = entry.getKey();
            final String value = entry.getValue();
            if (key.equals(DgeRecordTag.VERSION.name())) {
                header.setVersion(value);
            } else if (key.equals(DgeRecordTag.EXPRESSION_FORMAT.name())) {
                try {
                    header.setExpressionFormat(DgeHeader.ExpressionFormat.valueOf(value));
                } catch (IllegalArgumentException e) {
                    log.info("Unrecognized expression format '" + value + "' in " + inputName);
                }
            } else {
                log.info("Unrecognized field name '" + key + "' in @DGE record of " + inputName);
            }
        }
    }

    private DgeHeaderLibrary parseLibrary(final String line, final String inputName) {
        final LinkedHashMap<String, String> fields = parseLine(line, LIBRARY_RECORD_LABEL, inputName);
        final String uei = fields.get(LibraryRecordTag.UEI.name());
        if (uei == null) {
            log.warn("Ignoring DGE @LIBRARY line missing UEI tag in " + inputName + ";" + line);
            return null;
        }
        final DgeHeaderLibrary ret = new DgeHeaderLibrary(uei);
        for (final Map.Entry<String, String> entry: fields.entrySet()) {
            final String key = entry.getKey();
            final String value = entry.getValue();
            try {
                switch (LibraryRecordTag.valueOf(key)) {
                    case INPUT:         ret.setInput(new File(value));    break;
                    case INPUT_DGE:      ret.setInputDge(new File(value)); break;
                    case PREFIX:        ret.setPrefix(value);             break;
                    case REFERENCE:     ret.setReference(new File(value));break;
                    case UEI:           break; // set in ctor
                }
            } catch (IllegalArgumentException e) {
                ret.setTag(key, value);
            }

        }
        return ret;
    }

    private DgeHeaderCommand parseCommand(String line, String inputName) {
        final LinkedHashMap<String, String> fields = parseLine(line, COMMAND_RECORD_LABEL, inputName);
        final String commandLine = fields.get(CommandRecordTag.CL.name());
        if (commandLine == null) {
            log.warn("Ignoring DGE @COMMAND line missing CL tag in " + inputName + ";" + line);
            return null;
        }
        return new DgeHeaderCommand(commandLine);
    }



    private LinkedHashMap<String, String> parseLine(final String line, final String recordType, final String inputName) {
        final String[] fields = line.split(FIELD_SEPARATOR);
        if (!recordType.equals(fields[0])) {
            throw new IllegalStateException("Line did not have expected record type");
        }
        final LinkedHashMap<String, String> ret = new LinkedHashMap<>();
        for (int i = 1; i < fields.length; ++i) {
            final String[] kvPair = fields[i].split(KV_SEPARATOR, 2);
            if (kvPair.length != 2) {
                log.warn("Strange DGE header line seen in " + inputName + ";" + line);
                return null;
            } else {
                final String key = kvPair[0];
                final String value = kvPair[1];
                if (key.isEmpty() || value.isEmpty()) {
                    log.warn("Strange key-value pair in DGE header line seen in " + inputName + ";" + line);
                } else if (ret.containsKey(key)) {
                    log.warn("Key " + key + " appears more than once in DGE Header lineseen in " + inputName + ";" + line);
                } else {
                    ret.put(key, value);
                }
            }
        }
        return ret;
    }
}
