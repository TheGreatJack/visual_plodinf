#!/bin/bash

bam_file="$1"
out_path="$2"
minimum_ad="$3"
minimum_depth="$4"
maximum_depth="$5"
minimum_baseQ="$6"
minimum_mapQ="$7"

# Replace the following line with the command you want to execute on each file

prefix=$(basename $bam_file .bam)

echo "START------------------------------"
echo $bam_file
echo $out_path/$prefix\_$pop\_allele_prop_m$minimum_ad\_d$minimum_depth.tsv

bcftools mpileup -Q $minimum_baseQ -q $minimum_mapQ --skip-all-unset 3 -a "FORMAT/AD" -Ou --no-reference $bam_file | bcftools filter -Ou -e 'INFO/QS[1] == 1 | INFO/DP <= 10' - | bcftools query -f '%CHROM\t%POS\t%ALT\t%DP\tQS:%QS[\t%AD]]\n' | python allele_counts_parser.py -m $minimum_ad -d $minimum_depth -M $maximum_depth > $out_path/$prefix.tsv

echo "Finished processing file: $prefix"