#! /usr/bin/bash

##########################################
########### input/output params ##########
##########################################

VCF_DIR="/data/CCBR/projects/ccbr720/dragen_vcfs"
VCF_SUFFIX=".hard-filtered.vcf"   ### NOTE: vcf2maf does NOT take compressed VCFs, so will have to gunzip if it's a gzipped file
OUT_DIR="maf"
TUM_NORM_SEP="_vs_"

##########################################
############# vcf2maf params #############
##########################################

############# hg38 #############
SPECIES="homo_sapiens" ## "VEPSPECIES" variable in hg38.json
GENOME="GRCh38" ## "VEPBUILD" variable in hg38.json
GENOME_LOC="/data/CCBR_Pipeliner/db/PipeDB/lib/GRCh38.d1.vd1.fa"  ## "VEPFASTA" variable in hg38.json
FILTER_SITES="/data/CCBR_Pipeliner/db/PipeDB/lib/ExAC_nonTCGA.r0.3.1.sites.vep.GRCh38.vcf.gz"  ## "VEPFILTERVCF" variable in hg38.json

############# hg19 #############
#SPECIES="homo_sapiens" ## "VEPSPECIES" variable in hg19.json
#GENOME="GRCh37" ## "VEPBUILD" variable in hg19.json
#GENOME_LOC="/data/CCBR_Pipeliner/db/PipeDB/lib/hs37d5.fa"  ## "VEPFASTA" variable in hg19.json
#FILTER_SITES="/fdb/VEP/88/cache/ExAC.r0.3.sites.vep.vcf.gz"  ## "VEPFILTERVCF" variable in hg19.json

############# mm10 #############
#SPECIES="mus_musculus" ## "VEPSPECIES" variable in mm10.json
#GENOME="GRCm38" ## "VEPBUILD" variable in mm10.json
#GENOME_LOC="/data/CCBR_Pipeliner/db/PipeDB/lib/genome.fa"  ## "VEPFASTA" variable in mm10.json
#FILTER_SITES="/data/CCBR_Pipeliner/db/PipeDB/lib/mm10_allstrains_dbSNP142.vcf.gz"  ## "VEPFILTERVCF" variable in mm10.json


SWARM_FILE="vcf2maf.swarm"
if [ -e $SWARM_FILE ]; then
    rm $SWARM_FILE
    touch $SWARM_FILE
fi

if [ ! -e $OUT_DIR ]; then
    mkdir -p $OUT_DIR
fi

VCF_FILES=($(ls $VCF_DIR/*$TUM_NORM_SEP*$VCF_SUFFIX))

for file_name in "${VCF_FILES[@]}"; do    
    
    ### Example file name convention: "Pt01_Tumor_vs_Pt01_Normal.mutect2.vcf"
    ###    In this case, $TUM_NORM_SEP="_vs_", $VCF_SUFFIX=".mutect2.vcf", $norm_id="Pt01_Normal", $tumor_id="Pt01_Tumor"
    
    ### This assumes that the sample ID in the VCF header matches the file name
    ###    Better solution may be to use 'bcftools query -l' to list the samples in the VCF for the "vcf-tumor-id" and "vcf-normal-id" arguments to vcf2maf.pl
    ###    And use the extracted info from the file name as the "tumor-id" and "normal-id" arguments
    
    pairs=($(basename -s $VCF_SUFFIX $file_name | sed -e s/$TUM_NORM_SEP/\\n/))  ## Cleans and splits the filename to extract tumor ID
    tumor_id=${pairs[0]}
    norm_id=${pairs[1]}  ### Could add checking to see if this variable is empty, then do not include as argument
    norm_id_vcf=$(basename -s $VCF_SUFFIX $file_name)
    
    out_file="$OUT_DIR/$(basename -s $VCF_SUFFIX $file_name).maf"
    MY_CMD="date; module load VEP/97; module load vcf2maf/1.6.16; \
           vcf2maf.pl --input-vcf $file_name --output-maf $out_file \
           --vep-path \$VEP_HOME --vep-data \$VEPCACHEDIR \
            --species $SPECIES --ncbi-build $GENOME --ref-fasta $GENOME_LOC --filter-vcf $FILTER_SITES \
           --vcf-tumor-id $tumor_id --tumor-id $tumor_id \
           --vcf-normal-id $norm_id_vcf --normal-id $norm_id \   #### THIS SHOULD BE OMITTED FOR TUMOR-ONLY SAMPLES
           --vep-forks 2; date"
    echo $MY_CMD >> $SWARM_FILE
done

echo "Submitting swarm job..."
swarm_cmd="swarm -f $SWARM_FILE -g 16 -t 2 --logdir swarm_out/ --partition quick"
echo $swarm_cmd
#eval $swarm_cmd

echo "Done!"
