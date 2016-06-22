/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2016 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
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
    private enum DgeRecordTag {Version, ExpressionFormat}
    private enum LibraryRecordTag {Input, InputDge, Reference, Uei, Prefix}

    public void encode(final Writer writer, final DgeHeader header) {
        writeLine(writer, buildFirstLine(header));
        for (final String uei : new IterableAdapter<>(header.iterateLibraries())) {
            writeLine(writer, buildLibraryLine(header.getLibrary(uei)));
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
        dgeRecord.addFieldIfNotNull(DgeRecordTag.Version, header.getVersion());
        dgeRecord.addFieldIfNotNull(DgeRecordTag.ExpressionFormat, header.getExpressionFormat());
        return dgeRecord.build();
    }

    private String buildLibraryLine(final DgeHeaderLibrary library) {
        final OutputRecordBuilder dgeRecord = new OutputRecordBuilder(LIBRARY_RECORD_LABEL);
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.Input, library.getInput());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.InputDge, library.getInputDge());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.Reference, library.getReference());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.Uei, library.getUei());
        dgeRecord.addFieldIfNotNull(LibraryRecordTag.Prefix, library.getPrefix());
        for (final Map.Entry<String, String> entry: new IterableAdapter<>(library.iterateOtherTags())) {
            dgeRecord.addFieldIfNotNull(entry.getKey(), entry.getValue());
        }
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
        try {
            if (!reader.markSupported()) {
                throw new IllegalArgumentException("reader.markSupported == false");
            }
            final DgeHeader ret = new DgeHeader();
            ret.setExpressionFormat(null);
            ret.setVersion(null);
            boolean first = true;
            while (isHeaderLine(reader)) {
                final String line = reader.readLine();
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
            if (key.equals(DgeRecordTag.Version.name())) {
                header.setVersion(value);
            } else if (key.equals(DgeRecordTag.ExpressionFormat.name())) {
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
        final String uei = fields.get(LibraryRecordTag.Uei.name());
        if (uei == null) {
            log.warn("Ignoring DGE @LIBRARY line missing Uei tag in " + inputName + ";" + line);
            return null;
        }
        final DgeHeaderLibrary ret = new DgeHeaderLibrary(uei);
        for (final Map.Entry<String, String> entry: fields.entrySet()) {
            final String key = entry.getKey();
            final String value = entry.getValue();
            try {
                switch (LibraryRecordTag.valueOf(key)) {
                    case Input:         ret.setInput(new File(value));    break;
                    case InputDge:      ret.setInputDge(new File(value)); break;
                    case Prefix:        ret.setPrefix(value);             break;
                    case Reference:     ret.setReference(new File(value));break;
                    case Uei:           break; // set in ctor
                }
            } catch (IllegalArgumentException e) {
                ret.setTag(key, value);
            }

        }
        return ret;
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
