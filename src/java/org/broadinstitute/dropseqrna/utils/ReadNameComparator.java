package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;

import java.util.Comparator;

public class ReadNameComparator implements Comparator<SAMRecord>{

	public ReadNameComparator() {
    }

	@Override
    public int compare(final SAMRecord rec1, final SAMRecord rec2) {
        final String s1 = rec1.getReadName();
        final String s2 = rec2.getReadName();

        if (s1 != null) {
            if (s2 == null)
                return 1;
			else
				return s1.compareTo(s2);
        } else if (s2 != null)
			return -1;
		else
			return 0;
    }

}
