#!/usr/bin/env bash

set -euo pipefail

### Configuration ###

cram_dir="path_to_cram_directory"
bcf_dir="path_to_bcf_directory"
pileup_script="path_to_pileup.sh"
genome="path_to_genome.fa.gz"
norepeats_bed="path_to_directory/norepeats.bed"
post_merging_script="path_to_scripts/post_merging.awk"

mkdir -p "$bcf_dir"

### Variant Calling on mapped sequence files (.cram) ###

for cram_file in "$cram_dir"/*.cram; do
	echo "Processing: $cram_file …"
		base_name=$(basename "$cram_file" .cram)
		out_file="$bcf_dir/${base_name}.bcf"
		"$pileup_script" -f "$genome" -Ob "$cram_file" > "$out_file"
done

### Merging sample BCF files and inbred line BCF files ###

bcftools merge --no-index --merge both --threads 10 -Ov -o -"$bcf_dir"/*.bcf ore.bcf hel.bcf sam.bcf|\
"$post_merging_script"|\
bcftools view -Ob > mpileup.bcf

### Filtering the merged BCF file ###
#1. Filtering positions with high and low depth
#2. Removing repeat regions
#3. Keeping only SNPs
#4. Removing SNPs within 5 base pairs of INDELs


avg_depth=$(bcftools query -f '%DP\n' mpileup.bcf | awk 'BEGIN{a=0}{a+=$1}END{printf ("%d\n", a/NR)}')
flt_expr="TYPE = \"snp\" & INFO/DP > $((avg_depth / 2)) & INFO/DP < $((avg_depth * 2))"
bcftools filter -g 5 -T "$norepeats_bed" -Ou mpileup.bcf | bcftools filter -i "$flt_expr" -Ou - |bcftools norm -m- -Oz - > mpileup_flt.vcf.gz

### Generate merged marker files ###

bcftools query -L mpileup_flt.vcf.gz |\
grep -E '^Dmel_' > marker_samp.txt
bcftools view -S marker_samp.txt -Oz -o marker_flt.vcf.gz mpileup_flt.vcf.gz

### Filter for markers ###

# Preparing marker information files
	# Marker information files can be generated using Process_mark.R. Each inbred line will have its own separate file
	#You would need the VCF file containing all the inbred lines (marker_flt.vcf.gz)
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' ore_marker_file.txt|sort -u > ore_marker_ids.txt
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' sam_marker_file.txt|sort -u > sam_marker_ids.txt
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' hel_marker_file.txt|sort -u > hel_marker_ids.txt

#Prepare merged sample file
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o mpileup_flt_ids.vcf.gz mpileup_flt.vcf.gz
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o marker_flt_ids.vcf.gz marker_flt.vcf.gz
