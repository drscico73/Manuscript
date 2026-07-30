#!/bin/bash

### Variant Calling on mapped sequence files (.cram) ###
cram_dir="path_to_cram_directory"
bcf_dir="path_to_bcf_directory"
pileup_script="path_to_pileup.sh"
genome="path_to_genome.fa.gz"

for file in "$cram_dir"/*.cram; do
	echo "Processing $file …"
		base_name=$(basename "$file" .cram)
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

### Generate merged marker files ###

bcftools query -L mpileup_flt.vcf.gz |\
grep -E '^Dmel_' > marker_samp.txt
bcftools view -S marker_samp.txt -Oz -o marker_flt.vcf.gz mpileup_flt.vcf.gz

### Filter for markers ###

# Preparing marker files
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' ore_marker_file.txt|sort -u |ore_marker_ids.txt
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' sam_marker_file.txt|sort -u |sam_marker_ids.txt
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' hel_marker_file.txt|sort -u |hel_marker_ids.txt

#Prepare merged sample file
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o mpileup_flt_ids.vcf.gz mpileup_flt.vcf.gz
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o marker_flt_ids.vcf.gz marker_flt.vcf.gz

#Filtering

bcftools view -i 'ID=@ore_marker_ids.txt' -Oz -o ore_mpileup_flt.vcf.gz mpileup_flt_ids.vcf.gz
bcftools view -i 'ID=@hel_marker_ids.txt' -Oz -o hel_mpileup_flt.vcf.gz mpileup_flt_ids.vcf.gz
bcftools view -i 'ID=@sam_marker_ids.txt' -Oz -o sam_mpileup_flt.vcf.gz mpileup_flt_ids.vcf.gz
bcftools view -i 'ID=@ore_marker_ids.txt' -Oz -o ore_marker_flt.vcf.gz marker_flt_ids.vcf.gz
bcftools view -i 'ID=@hel_marker_ids.txt' -Oz -o hel_marker_flt.vcf.gz marker_flt_ids.vcf.gz
bcftools view -i 'ID=@sam_marker_ids.txt' -Oz -o sam_marker_flt.vcf.gz marker_flt_ids.vcf.gz

### Preparing sample-specific files ###

bcftools query -l mpileup_flt.vcf.gz | grep -E '^HO-|^Dmel_Helsinki_inbred_F11$|^Dmel_OregonR_non-inbred$' > ho_samples.txt
bcftools query -l mpileup_flt.vcf.gz | grep -E '^HS-|^Dmel_Helsinki_inbred_F11$|^Dmel_Samarkand_inbred$' > hs_samples.txt
bcftools query -l mpileup_flt.vcf.gz | grep -E '^SO-|^Dmel_Samarkand_inbred$|^Dmel_OregonR_non-inbred$' > so_samples.txt

bcftools view -S ho_samples.txt -Oz -o hel_ho_mpileup.vcf.gz hel_mpileup_flt.vcf.gz 
bcftools view -S hs_samples.txt -Oz -o hel_hs_mpileup.vcf.gz hel_mpileup_flt.vcf.gz
bcftools view -S so_samples.txt -Oz -o ore_so_mpileup.vcf.gz ore_mpileup_flt.vcf.gz
bcftools view -S ho_samples.txt -Oz -o ore_ho_mpileup.vcf.gz ore_mpileup_flt.vcf.gz
bcftools view -S hs_samples.txt -Oz -o sam_hs_mpileup.vcf.gz sam_mpileup_flt.vcf.gz
bcftools view -S so_samples.txt -Oz -o sam_so_mpileup.vcf.gz sam_mpileup_flt.vcf.gz

bcftools view -s ^Dmel_Samarkand_inbred -Oz -o hel_ho_markers.vcf hel_marker_flt.vcf.gz
bcftools view -s ^Dmel_OregonR_non-inbred -Oz -o hel_hs_markers.vcf hel_marker_flt.vcf.gz
bcftools view -s ^Dmel_OregonR_non-inbred -Oz -o sam_hs_markers.vcf sam_marker_flt.vcf.gz
bcftools view -s ^Dmel_Helsinki_inbred_F11 -Oz -o sam_so_markers.vcf sam_marker_flt.vcf.gz
bcftools view -s ^Dmel_Helsinki_inbred_F11 -Oz -o ore_so_markers.vcf ore_marker_flt.vcf.gz
bcftools view -s ^Dmel_Samarkand_inbred -Oz -o ore_ho_markers.vcf ore_marker_flt.vcf.gz
