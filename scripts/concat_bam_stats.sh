#!/bin/bash

# usage: bash concat_bam_stats.sh /path/to/stats_dir /path/to/output

STAT_DIR=$1
OUTPUT=$2

cd $STAT_DIR

# First just fetch the header
file1=$(ls *_collect_wgs_metrics.txt | head -n1)
grep "## METRIC" $file1 > $OUTPUT
grep "CATEGORY" $file1 >> $OUTPUT

for file in $(ls *_collect_wgs_metrics.txt)
do

  # Assign the sample name
  sample_name=$(echo $file | sed 's/_collect_wgs_metrics.txt//g')

  # And print to our output
  whole_genome_res=$(grep -w "WHOLE_GENOME" $file)
  non_zero_res=$(grep -w "NON_ZERO_REGIONS" $file)
  echo -e "${whole_genome_res}\t${sample_name}" >> $OUTPUT
  echo -e "${non_zero_res}\t${sample_name}" >> $OUTPUT
done

echo ">>> ALL FINISHED
"
