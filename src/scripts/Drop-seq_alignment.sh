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
ncores=1
bead_repair=0
keep_intermediates=0
pipeline=0


usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Drop-seq tagging, trimming and alignment.

-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-o <outputdir>      : Where to write output bam.  Default: '$outdir'.
-s <STAR_path>      : Full path of STAR.  Default: '$star_executable'.
-n <ncores>         : Number of cores to run.  Default: $ncores.
-b                  : Do bead repair.  Not needed for 10X libraries, but recommended for Drop-seq chemistry.
-e                  : Echo commands instead of executing them.
-k                  : Keep intermediate files.
-p                  : Run in pipeline mode.
-v                  : Run in verbose mode.
-h                  : Print usage and exit.
EOF
}


set -e

# Note: This script uses aliases to define commands containing quoted
# parameters for later use. By default, Bash only expands these in interactive
# shells. The following command enables that feature for the script, if
# necessary.
# shellcheck disable=SC3044 # `shopt` is only called for non-POSIX shells.
! command -v shopt > /dev/null || shopt -s expand_aliases

# Unset all variables to be set based on parameters that have not been
# initialized previously in this script (or a script sourced above) 'in order
# to ensure they're passed on command line rather than inherited from
# somewhere'.
# See https://github.com/broadinstitute/Drop-seq/pull/412#discussion_r1569231368
# for the corresponding discussion.
unset \
  genomedir \
  reference

while getopts ':g:r:o:s:n:bekpvh' options; do
  case $options in
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    o ) outdir=$OPTARG;;
    s ) star_executable=$OPTARG;;
    n ) ncores=$OPTARG;;
    b ) bead_repair=1;;
    e ) ECHO='echo';;
    k ) keep_intermediates=1;;
    p ) pipeline=1;;
    v ) verbose=1;;
    h ) usage
          exit 0;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $((OPTIND - 1))

check_set "$genomedir" 'Genome directory' '-g'
check_set "$reference" 'Reference fasta'  '-r'

check_TMPDIR

TMPDIR=$(mktemp -d)
echo "Writing temporary files to $TMPDIR" 

if [ "$#" -ne 1 ]
then error_exit 'Incorrect number of arguments'
fi

if [ "$star_executable" != 'STAR' ]
then if [ ! -x "$star_executable" ] || [ ! -f "$star_executable" ]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit 'STAR executable must be on the path'
fi

if [ "$pipeline" -ne 0 ]
then if [ -n "$ECHO" ]
     then error_exit 'Pipeline mode (-p flag) not supported in combination with echo mode (-e flag).'
     elif [ "$verbose" -ne 0 ]
     then error_exit 'Pipeline mode (-p flag) not supported in combination with verbose mode (-v flag).'
     fi
fi

reference_basename=$(basename "$(basename "$reference" .gz)" .fasta)
gene_intervals=$(dirname "$reference")/"$reference_basename".genes.intervals
refflat=$(dirname "$reference")/"$reference_basename".refFlat

unmapped_bam=$1
tagged_unmapped_bam=$TMPDIR/unaligned_mc_tagged_polyA_filtered.bam
aligned_sam=$TMPDIR/star.Aligned.out.sam
aligned_sorted_bam=$TMPDIR/aligned.sorted.bam


# Setup intermediate file cleanup

set --
# shellcheck disable=SC2142 # This is just for readability.
alias mark_file_as_intermediate='set -- "$@" '
# shellcheck disable=SC2317 # This function is called via a trap (see below).
cleanup_intermediates() {
    exit_code=$?
    if [ "$keep_intermediates" -eq 1 ]
    then
      [ "$verbose" -ne 1 ] || echo 'Keeping intermediate files.'
      exit "$exit_code"
    fi
    if [ "$exit_code" -ne 0 ]
    then
      [ "$verbose" -ne 1 ] || echo 'An error occurred...keeping intermediate files.'
      exit "$exit_code"
    fi
    [ "$verbose" -ne 1 ] || echo 'Finished successfully...removing intermediate files...'
    for intermediate_file do
      [ "$verbose" -ne 1 ] || echo "...removing intermediate file '$intermediate_file'."
      rm -- "$intermediate_file"
    done
}
trap 'cleanup_intermediates "$@"' EXIT


# Stage 1: pre-alignment tag and trim

# cellular tag
alias tag_cells='invoke_dropseq \
  TagBamWithReadSequenceExtended \
    BARCODED_READ=1 \
    BASE_QUALITY=10 \
    BASE_RANGE=1-12 \
    DISCARD_READ=false \
    INPUT="$unmapped_bam" \
    NUM_BASES_BELOW_QUALITY=1 \
    SUMMARY="$outdir"/unaligned_tagged_Cellular.bam_summary.txt \
    TAG_NAME=XC'

# molecular tag
alias tag_molecules='invoke_dropseq \
  TagBamWithReadSequenceExtended \
    BARCODED_READ=1 \
    BASE_QUALITY=10 \
    BASE_RANGE=13-20 \
    DISCARD_READ=true \
    NUM_BASES_BELOW_QUALITY=1 \
    SUMMARY="$outdir"/unaligned_tagged_Molecular.bam_summary.txt \
    TAG_NAME=XM'

# quality filter
alias filter_bam='invoke_dropseq \
  FilterBam \
    TAG_REJECT=XQ'

# read trimming
alias trim_start='invoke_dropseq \
  TrimStartingSequence \
    MISMATCHES=0 \
    NUM_BASES=5 \
    OUTPUT_SUMMARY="$outdir"/adapter_trimming_report.txt \
    SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG'

alias trim_polya='invoke_dropseq \
  PolyATrimmer \
    MISMATCHES=0 NUM_BASES=6 \
    NEW=true \
    OUTPUT="$tagged_unmapped_bam" \
    OUTPUT_SUMMARY="$outdir"/polyA_trimming_report.txt'

if [ "$pipeline" -ne 0 ]
then
  trap '{ exit 1; }' USR1
  pipe_fail () { kill -s USR1 "$$"; }
  { tag_cells \
      COMPRESSION_LEVEL=0 \
      OUTPUT=/dev/stdout \
    || pipe_fail; } \
  | { tag_molecules \
        COMPRESSION_LEVEL=0 \
        INPUT=/dev/stdin \
        OUTPUT=/dev/stdout \
    || pipe_fail; } \
  | { filter_bam \
        COMPRESSION_LEVEL=0 \
        INPUT=/dev/stdin \
        OUTPUT=/dev/stdout \
    || pipe_fail; } \
  | { trim_start \
        COMPRESSION_LEVEL=0 \
        INPUT=/dev/stdin \
        OUTPUT=/dev/stdout \
    || pipe_fail; } \
  | trim_polya \
      INPUT=/dev/stdin
else
  tag_cells \
    OUTPUT="$TMPDIR"/unaligned_tagged_Cell.bam
  mark_file_as_intermediate "$TMPDIR"/unaligned_tagged_Cell.bam

  tag_molecules \
    INPUT="$TMPDIR"/unaligned_tagged_Cell.bam \
    OUTPUT="$TMPDIR"/unaligned_tagged_CellMolecular.bam
  mark_file_as_intermediate "$TMPDIR"/unaligned_tagged_CellMolecular.bam

  filter_bam \
    INPUT="$TMPDIR"/unaligned_tagged_CellMolecular.bam \
    OUTPUT="$TMPDIR"/unaligned_tagged_filtered.bam
  mark_file_as_intermediate "$TMPDIR"/unaligned_tagged_filtered.bam

  trim_start \
    INPUT="$TMPDIR"/unaligned_tagged_filtered.bam \
    OUTPUT="$TMPDIR"/unaligned_tagged_trimmed_smart.bam
  mark_file_as_intermediate "$TMPDIR"/unaligned_tagged_trimmed_smart.bam

  trim_polya \
    INPUT="$TMPDIR"/unaligned_tagged_trimmed_smart.bam
fi
mark_file_as_intermediate "$tagged_unmapped_bam"


# Stage 2: alignment

alias sam_to_fastq='invoke_picard \
  SamToFastq \
    INPUT="$TMPDIR"/unaligned_mc_tagged_polyA_filtered.bam'

alias star_align='$ECHO \
  $star_executable \
    --genomeDir "$genomedir" \
    --outFileNamePrefix "$TMPDIR"/star. \
    --runThreadN $ncores'

if [ "$pipeline" -ne 0 ]
then
  { sam_to_fastq \
      COMPRESSION_LEVEL=0 \
      FASTQ=/dev/stdout \
    || pipe_fail; } \
  | star_align \
      --readFilesIn /dev/stdin
else
  sam_to_fastq \
    FASTQ="$TMPDIR"/unaligned_mc_tagged_polyA_filtered.fastq
  mark_file_as_intermediate "$TMPDIR"/unaligned_mc_tagged_polyA_filtered.fastq

  star_align \
    --readFilesIn "$TMPDIR"/unaligned_mc_tagged_polyA_filtered.fastq
fi
mark_file_as_intermediate "$aligned_sam"

$ECHO mv "$TMPDIR"/star.Log.final.out "$outdir"


# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)

alias sort_star='invoke_picard \
  SortSam \
    SORT_ORDER=queryname \
    INPUT="$aligned_sam" \
    OUTPUT="$aligned_sorted_bam" \
    TMP_DIR="$TMPDIR"'

sort_star
mark_file_as_intermediate "$aligned_sorted_bam"


# Stage 4: merge and tag aligned reads

alias merge_bam='invoke_picard \
  MergeBamAlignment \
    ALIGNED_BAM="$aligned_sorted_bam" \
    CLIP_ADAPTERS=false \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false \
    REFERENCE_SEQUENCE="$reference" \
    TMP_DIR="$TMPDIR" \
    UNMAPPED_BAM="$tagged_unmapped_bam"'

alias tag_with_interval='invoke_dropseq \
  TagReadWithInterval \
    INTERVALS="$gene_intervals" \
    TAG=XG \
    TMP_DIR="$TMPDIR"'

alias tag_with_gene='invoke_dropseq \
  TagReadWithGeneFunction \
    ANNOTATIONS_FILE="$refflat" \
    OUTPUT="$TMPDIR"/function_tagged.bam'

if [ "$pipeline" -ne 0 ]
then
  { merge_bam \
      COMPRESSION_LEVEL=0 \
      OUTPUT=/dev/stdout \
    || pipe_fail; } \
  | { tag_with_interval \
        COMPRESSION_LEVEL=0 \
        INPUT=/dev/stdin \
        OUTPUT=/dev/stdout \
      || pipe_fail; } \
  | tag_with_gene \
      INPUT=/dev/stdin
else
  merge_bam \
    OUTPUT="$TMPDIR"/merged.bam
  mark_file_as_intermediate "$TMPDIR"/merged.bam

  tag_with_interval \
    INPUT="$TMPDIR"/merged.bam \
    OUTPUT="$TMPDIR"/gene_tagged.bam
  mark_file_as_intermediate "$TMPDIR"/gene_tagged.bam

  tag_with_gene \
    INPUT="$TMPDIR"/gene_tagged.bam
fi


# Stage 5: bead repair

alias detect_subs_errors='invoke_dropseq \
  DetectBeadSubstitutionErrors \
    INPUT="$TMPDIR"/function_tagged.bam \
    MIN_UMIS_PER_CELL=20 \
    OUTPUT_REPORT="$outdir"/substitution_error_report.txt \
    TMP_DIR="$TMPDIR"'

alias detect_synthesis_errors='invoke_dropseq \
  DetectBeadSynthesisErrors \
    CREATE_INDEX=true \
    MIN_UMIS_PER_CELL=20 \
    OUTPUT="$outdir"/final.bam \
    OUTPUT_STATS="$outdir"/synthesis_error_stats.txt \
    REPORT="$outdir"/synthesis_error_report.txt \
    SUMMARY="$outdir"/synthesis_error_summary.txt \
    TMP_DIR="$TMPDIR"'

if [ "$bead_repair" -ne 0 ]
then
  mark_file_as_intermediate "$TMPDIR"/function_tagged.bam

  # Note: For some reason, piping the output of the fist command to the second
  # results in an empty BAM file. Thus, at least for the time being, pipeline
  # mode does not affect this stage.
  detect_subs_errors \
    OUTPUT="$TMPDIR"/substitution_repaired.bam
  mark_file_as_intermediate "$TMPDIR"/substitution_repaired.bam

  detect_synthesis_errors \
    INPUT="$TMPDIR"/substitution_repaired.bam
else
  $ECHO mv "$TMPDIR"/function_tagged.bam "$outdir"/final.bam
fi


echo 'Completed successfully.'
