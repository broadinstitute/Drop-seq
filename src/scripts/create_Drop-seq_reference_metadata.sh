#!/bin/sh

# MIT License
#
# Copyright 2018 Broad Institute
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


. "$(dirname "$0")"/defs.sh

progname=$(basename "$0")

outdir=$(pwd)
star_executable=$(which STAR 2> /dev/null)
bgzip_executable=$(which bgzip 2> /dev/null)
samtools_executable=$(which samtools 2> /dev/null)


usage () {
    cat >&2 <<EOF
USAGE: $progname [options]
Create Drop-seq reference metadata bundle.

-n <name>                  : Name for reference metadata set to be created.  Required.
-r <referencefasta>        : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-s <species>               : Species.  Required.
-g <gtf>                   : Gene annotation file.  Required.
-f <filtered-gene-biotype> : Annotations with the given gene biotype will be filtered.  Multiple
                             values may be specified by using this argument more than once, and/or
                             by providing a comma-separated list.  Use ValidateReference command to
                             see the gene biotypes in your GTF in order to decide what to exclude.
                             Default: No gene biotypes are filtered.
-o <outputdir>             : Where to write output bam.  Default: '$outdir'.
-a <STAR_path>             : Full path of STAR.  Default: '$star_executable'.
-b <bgzip_path>            : Full path of bgzip.  Default: '$bgzip_executable'.
-i <samtools_path>         : Full path of samtools.  Default: '$samtools_executable'.
-v                         : Run in verbose mode.
-e                         : Merely echo commands instead of executing.
-h                         : Print usage and exit.
EOF
}


set -e

# Unset all variables to be set based on parameters that have not been
# initialized previously in this script (or a script sourced above) 'in order
# to ensure they're passed on command line rather than inherited from
# somewhere'.
# See https://github.com/broadinstitute/Drop-seq/pull/412#discussion_r1569231368
# for the corresponding discussion.
unset \
  reference_name \
  reference_fasta \
  species \
  gtf \
  filtered_gene_biotypes

while getopts ':n:r:s:g:f:o:a:b:i:veh' options; do
  case $options in
    n ) reference_name=$OPTARG;;
    r ) reference_fasta=$OPTARG;;
    s ) species=$OPTARG;;
    g ) gtf=$OPTARG;;
    f ) filtered_gene_biotypes="${filtered_gene_biotypes+$filtered_gene_biotypes }G=$OPTARG";;
        # support repeating -f: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    o ) outdir=$OPTARG;;
    a ) star_executable=$OPTARG;;
    b ) bgzip_executable=$OPTARG;;
    i ) samtools_executable=$OPTARG;;
    v ) verbose=1;;
    e ) ECHO='echo';;
    h ) usage
          exit 0;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $((OPTIND - 1))

check_TMPDIR

check_set "$reference_name" 'Reference name' '-n'
check_set "$reference_fasta" 'Reference fasta' '-r'
check_set "$species" 'Species' '-s'
check_set "$gtf" 'Gene annotation file' '-g'
check_set "$star_executable" 'STAR path' '-s' 'if STAR is not on PATH'
check_set "$samtools_executable" 'samtools path' '-i' 'if samtools is not on PATH'
check_set "$bgzip_executable" 'bgzip path' '-b' 'if bgzip is not on PATH'


output_fasta=$outdir/$reference_name.fasta
sequence_dictionary=$outdir/$reference_name.dict
output_gtf=$outdir/$reference_name.gtf
reduced_gtf=$outdir/$reference_name.reduced.gtf


invoke_picard NormalizeFasta INPUT="$reference_fasta" OUTPUT="$output_fasta"
if [ -e "$sequence_dictionary" ]
then $ECHO rm "$sequence_dictionary"
fi
invoke_picard CreateSequenceDictionary REFERENCE="$output_fasta" OUTPUT="$sequence_dictionary" SPECIES="$species"
invoke_dropseq FilterGtf GTF="$gtf" SEQUENCE_DICTIONARY="$sequence_dictionary" OUTPUT="$output_gtf" $filtered_gene_biotypes
invoke_dropseq ConvertToRefFlat ANNOTATIONS_FILE="$output_gtf" SEQUENCE_DICTIONARY="$sequence_dictionary" OUTPUT="$outdir/$reference_name".refFlat
invoke_dropseq ReduceGtf GTF="$output_gtf" SEQUENCE_DICTIONARY="$sequence_dictionary" OUTPUT="$reduced_gtf"
invoke_dropseq CreateIntervalsFiles SEQUENCE_DICTIONARY="$sequence_dictionary" REDUCED_GTF="$reduced_gtf" PREFIX="$reference_name" \
           OUTPUT="$outdir"
$ECHO mkdir -p "$outdir"/STAR
check_invoke "$star_executable" --runMode genomeGenerate --genomeDir "$outdir"/STAR --genomeFastaFiles "$output_fasta" \
           --sjdbGTFfile "$output_gtf" --sjdbOverhang 100
check_invoke "$bgzip_executable" "$output_fasta"
check_invoke "$samtools_executable" faidx "$output_fasta".gz


echo 'Reference metadata created sucessfully.'
