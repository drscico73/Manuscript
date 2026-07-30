#!/usr/bin/env bash

set -euo pipefail

### Configuration ###
output_dir="path_to_results_directory"
mkdir -p "$output_dir"

# Main output directory 
output_dir="path_to_results_directory" 

# Input files 
ore_marker_file="path_to_directory/ore_marker_file.txt" 
sam_marker_file="path_to_directory/sam_marker_file.txt" 
hel_marker_file="path_to_directory/hel_marker_file.txt" 
mpileup_flt="path_to_directory/mpileup_flt.vcf.gz" 
marker_flt="path_to_directory/marker_flt.vcf.gz"

# Output directories
marker_ids_dir="$output_dir/marker_ids" 
ids_dir="$output_dir/with_ids" 
filtered_dir="$output_dir/filtered" 
samples_dir="$output_dir/sample_specific"
mkdir -p "$marker_ids_dir" "$ids_dir" "$filtered_dir" "$samples_dir"

###################################################
### Files needed ### 
# 1. OregonR marker information file: ore_marker_file.txt 
# 2. Samarkand marker information file: sam_marker_file.txt 
# 3. Helsinki marker information file: hel_marker_file.txt 
# 4. Merged and filtered sample file: mpileup_flt.vcf.gz 
# 5. Merged and filtered marker file: marker_flt.vcf.gz 
# 
# Marker information files can be generated using Process_mark.R. 
# Each inbred line should have its own separate marker information file.
###################################################

### Filter for markers ###

# Preparing marker information files
	# Marker information files can be generated using Process_mark.R. Each inbred line will have its own separate file
	#You would need the VCF file containing all the inbred lines (marker_flt.vcf.gz)

echo "Preparing marker ID files..."

awk 'NR > 1 {printf $1":"$2":"$3":"$4}' "$ore_marker_file"|sort -u > "$marker_ids_dir/ore_marker_ids.txt"
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' "$sam_marker_file"|sort -u > "$marker_ids_dir/sam_marker_ids.txt"
awk 'NR > 1 {printf $1":"$2":"$3":"$4}' "$hel_marker_file"|sort -u > "$marker_ids_dir/hel_marker_ids.txt"

# Adding IDs to merged sample and marker VCF files

echo "Adding variant IDs to VCF files..."

bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$ids_dir/mpileup_flt_ids.vcf.gz" "$mpileup_flt"
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$ids_dir/marker_flt_ids.vcf.gz" "$marker_flt"

# Filter merged sample and marker files using marker IDs

echo "Filtering sample and marker VCF by marker IDs..."

bcftools view -i 'ID=@${marker_ids_dir}/ore_marker_ids.txt' -Oz -o "$filtered_dir/ore_mpileup_flt.vcf.gz" "$ids_dir/mpileup_flt_ids.vcf.gz"
bcftools view -i 'ID=@${marker_ids_dir}/hel_marker_ids.txt' -Oz -o "$filtered_dir/hel_mpileup_flt.vcf.gz" "$ids_dir/mpileup_flt_ids.vcf.gz"
bcftools view -i 'ID=@${marker_ids_dir}/sam_marker_ids.txt' -Oz -o "$filtered_dir/sam_mpileup_flt.vcf.gz" "$ids_dir/mpileup_flt_ids.vcf.gz"
bcftools view -i 'ID=@${marker_ids_dir}/ore_marker_ids.txt' -Oz -o "$filtered_dir/ore_marker_flt.vcf.gz" "$ids_dir/marker_flt_ids.vcf.gz"
bcftools view -i 'ID=@${marker_ids_dir}/hel_marker_ids.txt' -Oz -o "$filtered_dir/hel_marker_flt.vcf.gz" "$ids_dir/marker_flt_ids.vcf.gz"
bcftools view -i 'ID=@${marker_ids_dir}/sam_marker_ids.txt' -Oz -o "$filtered_dir/sam_marker_flt.vcf.gz" "$ids_dir/marker_flt_ids.vcf.gz"

### Preparing sample-specific files ###

echo "Preparing sample lists..."

bcftools query -l "$mpileup_flt" | grep -E '^HO-|^Dmel_Helsinki_inbred_F11$|^Dmel_OregonR_non-inbred$' > "$samples_dir/ho_samples.txt"
bcftools query -l "$mpileup_flt" | grep -E '^HS-|^Dmel_Helsinki_inbred_F11$|^Dmel_Samarkand_inbred$' > "$samples_dir/hs_samples.txt"
bcftools query -l "$mpileup_flt" | grep -E '^SO-|^Dmel_Samarkand_inbred$|^Dmel_OregonR_non-inbred$' > "$samples_dir/so_samples.txt"

echo "Generating sample-specific mpileup files..."

bcftools view -S "$samples_dir/ho_samples.txt" -Oz -o "$samples_dir/hel_ho_mpileup.vcf.gz" "$filtered_dir/hel_mpileup_flt.vcf.gz"
bcftools view -S "$samples_dir/hs_samples.txt" -Oz -o "$samples_dir/hel_hs_mpileup.vcf.gz" "$filtered_dir/hel_mpileup_flt.vcf.gz"
bcftools view -S "$samples_dir/so_samples.txt" -Oz -o "$samples_dir/ore_so_mpileup.vcf.gz" "$filtered_dir/ore_mpileup_flt.vcf.gz"
bcftools view -S "$samples_dir/ho_samples.txt" -Oz -o "$samples_dir/ore_ho_mpileup.vcf.gz" "$filtered_dir/ore_mpileup_flt.vcf.gz"
bcftools view -S "$samples_dir/hs_samples.txt" -Oz -o "$samples_dir/sam_hs_mpileup.vcf.gz" "$filtered_dir/sam_mpileup_flt.vcf.gz"
bcftools view -S "$samples_dir/so_samples.txt" -Oz -o "$samples_dir/sam_so_mpileup.vcf.gz" "$filtered_dir/sam_mpileup_flt.vcf.gz"

echo "Generating sample-specific marker files..."

bcftools view -s ^Dmel_Samarkand_inbred -Oz -o "$samples_dir/hel_ho_markers.vcf.gz" "$filtered_dir/hel_marker_flt.vcf.gz"
bcftools view -s ^Dmel_OregonR_non-inbred -Oz -o "$samples_dir/hel_hs_markers.vcf.gz" "$filtered_dir/hel_marker_flt.vcf.gz"
bcftools view -s ^Dmel_OregonR_non-inbred -Oz -o "$samples_dir/sam_hs_markers.vcf.gz" "$filtered_dir/sam_marker_flt.vcf.gz"
bcftools view -s ^Dmel_Helsinki_inbred_F11 -Oz -o "$samples_dir/sam_so_markers.vcf.gz" "$filtered_dir/sam_marker_flt.vcf.gz"
bcftools view -s ^Dmel_Helsinki_inbred_F11 -Oz -o "$samples_dir/ore_so_markers.vcf.gz" "$filtered_dir/ore_marker_flt.vcf.gz"
bcftools view -s ^Dmel_Samarkand_inbred -Oz -o "$samples_dir/ore_ho_markers.vcf.gz" "$filtered_dir/ore_marker_flt.vcf.gz"

###done

echo "Processing completed"


