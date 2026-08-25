#!/usr/bin/env bash

set -euo pipefail

### Configuration ###

cram_dir="path_to_cram_directory"
bcf_dir="path_to_bcf_directory"
pileup_script="path_to_pileup.sh"
genome="path_to_genome.fa.gz"
norepeats_bed="path_to_directory/norepeats.bed"
post_merging_script="path_to_scripts/post_merging.awk"
output_dir="path_to_results_directory" 

mkdir -p "$bcf_dir"
mkdir -p "$output_dir"

### Variant Calling on mapped sequence files (.cram) ###

for cram_file in "$cram_dir"/*.cram; do
	echo "Processing: $cram_file …"
		base_name=$(basename "$cram_file" .cram)
		out_file="$bcf_dir/${base_name}.bcf"
		"$pileup_script" -f "$genome" -Ob "$cram_file" > "$out_file"
done

### Merging sample BCF files and inbred line BCF files ###

echo "Merging sample files"

bcftools merge --no-index --merge both --threads 10 -Ov -o -"$bcf_dir"/*.bcf ore.bcf hel.bcf sam.bcf|\
"$post_merging_script"|\
bcftools view -Ob > "$output_dir/mpileup.bcf"

### Filtering the merged BCF file ###
#1. Filtering positions with high and low depth
#2. Removing repeat regions
#3. Keeping only SNPs
#4. Removing SNPs within 5 base pairs of INDELs


avg_depth=$(bcftools query -f '%DP\n' "$output_dir/mpileup.bcf" | awk 'BEGIN{a=0}{a+=$1}END{printf ("%d\n", a/NR)}')
flt_expr="TYPE = \"snp\" & INFO/DP > $((avg_depth / 2)) & INFO/DP < $((avg_depth * 2))"

echo "Average depth: "$avg_depth""
bcftools filter -g 5 -T "$norepeats_bed" -Ou mpileup.bcf | bcftools filter -i "$flt_expr" -Ou - |bcftools norm -m- -Oz - > "$output_dir/mpileup_flt.vcf.gz"

### Generate merged marker files ###

echo "Generating marker file" 

bcftools query -L "$output_dir/mpileup_flt.vcf.gz" |\
grep -E '^Dmel_' > "$output_dir/marker_samp.txt"
bcftools view -S marker_samp.txt -Oz -o "$output_dir/marker_flt.vcf.gz" "$output_dir/mpileup_flt.vcf.gz"

####Done
echo "Processing completed."
