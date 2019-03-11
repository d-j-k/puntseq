#!/usr/bin/env sh

is_multiplexed="$1"
fastq="$2"
out="$3"
output="$4"
threads="$5"
check_reads="$6"
output_fmt="$7"

if [[ $is_multiplexed = "True" ]]; then
    porechop --input "$fastq" \
      --barcode_dir "$out" \
      --threads "$threads" \
      --check_reads "$check_reads" \
      --discard_middle \
      --discard_unassigned \
      --format "$output_fmt"

    echo "Porechop demultiplex and adapter trimming finished."
else
    porechop --input "$fastq" \
      --output "$out" \
      --threads "$threads" \
      --check_reads "$check_reads" \
      --discard_middle \
      --format "$output_fmt"

      echo "Porechop adapter trimming finished."
fi

touch "$output"
