#!/bin/bash

### Variant Calling on mapped sequence files (.cram) ###
cram_dir="path_to_cram_directory"
bcf_dir="path_to_bcf_directory"
pileup_script="path_to_pileup.sh"
genome="path_to_genome.fa.gz"

for file in "$cram_dir"/*.cram; do
	echo "Processing $file …"
		base_name=$(basename “$file” .cram)
		out_file="$bcf_dir/${base_name}.bcf"
		"$pileup_script" -f "$genome" -Ob "$file" > "$out_file"
done

### Merging Sample files and inbred line bcf files ###

bcftools merge --no-index --merge both --threads 10 -Ov -o -"$bcf_dir"/*.bcf ore.bcf hel.bcf sam.bcf|\
path_to_scripts/post_merging.awk|\
bcftools view -Ob > mpileup.bcf

### Filtering the merged bcf file ###
#1. Filtering positions with high and low depth
#2. Removing repeat regions
#3. Keeping only SNPs
#4. Removing SNPs within 5 base pairs of INDELs

norepeats_bed="path_to_directory/norepeats.bed"
avg_depth=$(bcftools query -f '%DP\n' mpileup.bcf | awk 'BEGIN{a=0}{a+=$1}END{printf ("%d\n", a/NR)}')
flt_expr="TYPE = \"snp\" & INFO/DP > $((avg_depth / 2)) & INFO/DP < $((avg_depth * 2))"
bcftools filter -g 5 -T "$norepeats_bed" -Ou mpileup.bcf | bcftools filter -i "$flt_expr" -Ou - |bcftools norm -m- -Oz - > mpileup_flt.vcf.gz

### Filter for markers ###

# Preparing marker files
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' ore_marker_file.txt|sort -u |ore_marker_ids.txt
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' sam_marker_file.txt|sort -u |sam_marker_ids.txt
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' hel_marker_file.txt|sort -u |hel_marker_ids.txt

#Prepare merged sample file
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o mpileup_flt_ids.vcf.gz mpileup_flt.vcf.gz

#
