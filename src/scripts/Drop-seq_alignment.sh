#!/usr/bin/env bash
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

source "$(dirname "$0")"/defs.sh

outdir=$(pwd)
genomedir=
reference=
star_executable=STAR
keep_intermediates=0
bead_repair=0
progname=$(basename "$0")

usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Drop-seq tagging, trimming and alignment

-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-b                  : Do bead repair.  Not needed for 10X libraries, but recommended for Drop-seq chemistry.  Default: disabled.
-e                  : Echo commands instead of executing them.
-k                  : Keep intermediate files
-v                  : verbose
EOF
}

set -e
# Fail if any of the commands in a pipeline fails
set -o pipefail

while getopts ":o:g:r:es:kvbh" options; do
  case $options in
    o ) outdir=$OPTARG;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    s ) star_executable=$OPTARG;;
    e ) ECHO="echo";;
    b ) bead_repair=1;;
    k ) keep_intermediates=1;;
    v ) verbose=1;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $((OPTIND - 1))

check_set "$genomedir" "Genome directory" "-g"
check_set "$reference" "Reference fasta"  "-r"

check_TMPDIR

TMPDIR=$(mktemp -d)
echo "Writing temporary files to $TMPDIR" 

if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi


if [[ "$star_executable" != "STAR" ]]
then if [[ ! ( -x $star_executable && -f $star_executable ) ]]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit "STAR executable must be on the path"
fi

reference_basename=$(basename "$(basename "$reference" .gz)" .fasta)
gene_intervals=$(dirname "$reference")/"$reference_basename".genes.intervals
refflat=$(dirname "$reference")/"$reference_basename".refFlat

unmapped_bam=$1
tagged_unmapped_bam=${TMPDIR}/unaligned_mc_tagged_polyA_filtered.bam
aligned_sam=${TMPDIR}/star.Aligned.out.sam
aligned_sorted_bam=${TMPDIR}/aligned.sorted.bam
files_to_delete=()


# Stage 1: pre-alignment tag and trim

# cellular tag
invoke_dropseq TagBamWithReadSequenceExtended SUMMARY="${outdir}"/unaligned_tagged_Cellular.bam_summary.txt \
  BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 \
  INPUT="${unmapped_bam}" OUTPUT="$TMPDIR"/unaligned_tagged_Cell.bam
files_to_delete+=("$TMPDIR/unaligned_tagged_Cell.bam")

# molecular tag
invoke_dropseq TagBamWithReadSequenceExtended SUMMARY="${outdir}"/unaligned_tagged_Molecular.bam_summary.txt \
  BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 \
  INPUT="$TMPDIR"/unaligned_tagged_Cell.bam OUTPUT="$TMPDIR"/unaligned_tagged_CellMolecular.bam
files_to_delete+=("$TMPDIR/unaligned_tagged_CellMolecular.bam")

# quality filter
invoke_dropseq FilterBam TAG_REJECT=XQ INPUT="$TMPDIR"/unaligned_tagged_CellMolecular.bam \
  OUTPUT="$TMPDIR"/unaligned_tagged_filtered.bam
files_to_delete+=("$TMPDIR/unaligned_tagged_filtered.bam")

# read trimming
invoke_dropseq TrimStartingSequence OUTPUT_SUMMARY="${outdir}"/adapter_trimming_report.txt \
  SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 INPUT="$TMPDIR"/unaligned_tagged_filtered.bam \
  OUTPUT="$TMPDIR"/unaligned_tagged_trimmed_smart.bam
files_to_delete+=("$TMPDIR/unaligned_tagged_trimmed_smart.bam")

invoke_dropseq PolyATrimmer OUTPUT="${tagged_unmapped_bam}" OUTPUT_SUMMARY="${outdir}"/polyA_trimming_report.txt \
  MISMATCHES=0 NUM_BASES=6 NEW=true INPUT="$TMPDIR"/unaligned_tagged_trimmed_smart.bam
files_to_delete+=("${tagged_unmapped_bam}")


# Stage 2: alignment
invoke_picard SamToFastq INPUT="${TMPDIR}"/unaligned_mc_tagged_polyA_filtered.bam \
  FASTQ="$TMPDIR"/unaligned_mc_tagged_polyA_filtered.fastq
files_to_delete+=("$TMPDIR/unaligned_mc_tagged_polyA_filtered.fastq")

$ECHO "$star_executable" --genomeDir "${genomedir}" --outFileNamePrefix "${TMPDIR}"/star. \
  --readFilesIn "$TMPDIR"/unaligned_mc_tagged_polyA_filtered.fastq
files_to_delete+=("${aligned_sam}")

# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
invoke_picard \
  SortSam INPUT="${aligned_sam}" OUTPUT="${aligned_sorted_bam}" SORT_ORDER=queryname TMP_DIR="${TMPDIR}"
files_to_delete+=("${aligned_sorted_bam}")

# Stage 4: merge and tag aligned reads
invoke_picard MergeBamAlignment REFERENCE_SEQUENCE="${reference}" UNMAPPED_BAM="${tagged_unmapped_bam}" \
  ALIGNED_BAM="${aligned_sorted_bam}" INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false CLIP_ADAPTERS=false \
  TMP_DIR="${TMPDIR}" OUTPUT="$TMPDIR"/merged.bam
files_to_delete+=("$TMPDIR/merged.bam")

invoke_dropseq TagReadWithInterval I="$TMPDIR"/merged.bam O="$TMPDIR"/gene_tagged.bam TMP_DIR="${TMPDIR}" \
  INTERVALS="${gene_intervals}" TAG=XG
files_to_delete+=("$TMPDIR/gene_tagged.bam")

invoke_dropseq TagReadWithGeneFunction O="${TMPDIR}"/function_tagged.bam ANNOTATIONS_FILE="${refflat}" \
  INPUT="$TMPDIR"/gene_tagged.bam

if [ "$bead_repair" -ne 0 ]
then
  files_to_delete+=("$TMPDIR/function_tagged.bam")

  # Stage 5: bead repair
  invoke_dropseq DetectBeadSubstitutionErrors INPUT="${TMPDIR}"/function_tagged.bam OUTPUT="${TMPDIR}"/substitution_repaired.bam \
    TMP_DIR="$TMPDIR" MIN_UMIS_PER_CELL=20 OUTPUT_REPORT="${outdir}"/substitution_error_report.txt
  files_to_delete+=("${TMPDIR}/substitution_repaired.bam")

  invoke_dropseq DetectBeadSynthesisErrors INPUT="${TMPDIR}"/substitution_repaired.bam MIN_UMIS_PER_CELL=20 \
    OUTPUT_STATS="${outdir}"/synthesis_error_stats.txt SUMMARY="${outdir}"/synthesis_error_summary.txt \
    REPORT="${outdir}"/synthesis_error_report.txt CREATE_INDEX=true TMP_DIR="$TMPDIR" OUTPUT="$outdir"/final.bam
else
  $ECHO mv "$TMPDIR"/function_tagged.bam "$outdir"/final.bam
fi

if [ "$keep_intermediates" -eq 0 ]
then $ECHO rm "${files_to_delete[@]}"
fi

echo "Completed successfully."

