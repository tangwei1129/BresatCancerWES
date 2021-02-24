vep \
 -i Sample_15120+Sample_15121.FINALmutect2.vcf \
 --vcf \
 -o example.vcf \
 --offline \
 --cache \
 --force_overwrite \
 --fork 16 \
 --dir_cache /data/tangw3/VEP/customDB \
 --species human \
 --assembly GRCh38 \
 --fasta  customDB/GRCh38.d1.vd1.fa \
 --af_gnomad \
 --custom  /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,variant_type,segdup

vcf2maf.pl \
  --input-vcf example.vcf \
  --output-maf example.maf \
  --inhibit-vep \
  --ref-fasta  customDB/GRCh38.d1.vd1.fa \
  --vep-forks 16



vcf2maf.pl --input-vcf $file_name --output-maf $out_file \
           --vep-path \$VEP_HOME --vep-data \$VEPCACHEDIR \
            --species $SPECIES --ncbi-build $GENOME --ref-fasta $GENOME_LOC --filter-vcf $FILTER_SITES \
           --vcf-tumor-id $tumor_id --tumor-id $tumor_id \
           --vcf-normal-id $norm_id_vcf --normal-id $norm_id \   #### THIS SHOULD BE OMITTED FOR TUMOR-ONLY SAMPLES
           --vep-forks 2;

### vcf2maf with VEP ####
vcf2maf.pl \
  --input-vcf $VEP_EXAMPLES/homo_sapiens_GRCh38.vcf \
  --output-maf example_38.maf \
  --vep-path $VEP_HOME \
  --vep-data $VEP_CACHEDIR \
  --ref-fasta $VEP_CACHEDIR/GRCh38.fa \
  --vep-forks 16 \
  --tumor-id Sample_15120 \
  --normal-id Sample_15121 \
  --ncbi-build GRCh38 





vcf2maf.pl \
  --input-vcf Sample_15120+Sample_15121.FINALmutect2.vcf \
  --output-maf test_37.maf \
  --vep-path $VEP_HOME \
  --vep-data $VEP_CACHEDIR \
  --vep-forks 16 \
  --tumor-id Sample_15120 \
  --normal-id Sample_15121 \
  --ncbi-build GRCh38 


 # --ref-fasta $VEP_CACHEDIR/GRCh38.fa \





vep  -i Sample_15120+Sample_15121.FINALmutect2.vcf -o example.vcf  --offline  --cache  --force_overwrite  --fork 16  --dir_cache /data/tangw3/VEP/customDB  --species human  --assembly GRCh38  --fasta  customDB/GRCh38.d1.vd1.fa  --everything --custom /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,variant_type,segdup  --vcf  --field "gnomADg_variant_type,gnomADg_AF_AFR"

species human
assembly GRCh38
no_progress
no_stats
sift b
ccds
uniprot
hgvs
symbol
numbers
domains
gene_phenotype
canonical
protein
biotype
uniprot
tsl
variant_class
shift_hgvs 1
check_existing
total_length
allele_number
no_escape
xref_refseq
failed 1
vcf
flag_pick_allele
pick_order canonical,tsl,biotype,rank,ccds,length
format vcf
af_1kg
af_esp
af_gnomad
vep-forks 16

cache
dir_cache /data/tangw3/VEP/customDB 
offline
fasta  customDB/GRCh38.d1.vd1.fa
show_cache_info
offline
force_overwrite



vcf2maf.pl   --input-vcf example.vcf   --output-maf example.maf   --inhibit-vep   --ref-fasta  customDB/GRCh38.d1.vd1.fa --tumor-id Sample_15120 --normal-id Sample_15121  --vep-forks 16


vep  -i Sample_15120+Sample_15121.FINALmutect2.vcf  --vcf  -o ex.vcf  --offline  --cache  --force_overwrite  --fork 16  --dir_cache /data/tangw3/VEP/customDB  --species human  --assembly GRCh38  --fasta  customDB/GRCh38.d1.vd1.fa --custom  /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomAD,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,variant_type,segdup


vcf2maf.pl   --input-vcf ex.vcf   --output-maf ex.maf   --inhibit-vep   --ref-fasta  customDB/GRCh38.d1.vd1.fa --tumor-id Sample_15120 --normal-id Sample_15121  --vep-forks 16


vep  -i Sample_15120+Sample_15121.FINALmutect2.vcf -o example.vcf \
--species human \
--assembly GRCh38 \
--no_progress \
--no_stats \
--sift b \
--ccds \
--uniprot \
--hgvs \
--symbol \
--numbers \
--domains \
--gene_phenotype \
--canonical \
--protein \
--biotype \
--uniprot \
--tsl \
--variant_class \
--shift_hgvs 1 \
--check_existing \
--total_length \
--allele_number \
--no_escape \
--xref_refseq \
--failed 1 \
--vcf \
--flag_pick_allele \
--pick_order canonical,tsl,biotype,rank,ccds,length \
--format vcf \
--af_1kg \
--af_esp \
--af_gnomad \
--fork 16 \
--cache \
--dir_cache /data/tangw3/VEP/customDB \
--offline \
--fasta  customDB/GRCh38.d1.vd1.fa \
--force_overwrite

vcf2maf.pl   --input-vcf example.vcf   --output-maf example.maf   --inhibit-vep   --ref-fasta  customDB/GRCh38.d1.vd1.fa --tumor-id Sample_15120 --normal-id Sample_15121  --vep-forks 16





vep --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir /fdb/VEP/101/cache --fasta customDB/GRCh38.d1.vd1.fa --format vcf --input_file Sample_15120+Sample_15121.FINALmutect2.vcf --output_file ./Sample_15120+Sample_15121.FINALmutect2.vep.vcf --offline --pubmed --fork 16 --polyphen b --af --af_1kg --af_esp --af_gnomad --regulatory --force_overwrite --pick --plugin ExAC,$VEP_CACHEDIR/ExAC.r0.3.1.sites.vep.vcf.gz --custom /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,variant_type,segdup







vep --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir /fdb/VEP/101/cache --fasta customDB/GRCh38.d1.vd1.fa --format vcf --input_file Sample_15120+Sample_15121.FINALmutect2.vcf --output_file ./Sample_15120+Sample_15121.FINALmutect2.vep.vcf --offline --pubmed --fork 16 --polyphen b --af --af_1kg --af_esp --af_gnomad --regulatory --force_overwrite --pick --plugin ExAC,/data/tangw3/VEP/customDB/ExAC.0.3.GRCh38.vcf.gz

--custom $VEP_CACHEDIR/gnomad.exomes.r2.0.2.sites.GRCh38.noVEP.vcf.gz,gnomADe,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS,segdup,lcr,variant_type
--custom $VEP_CACHEDIR/gnomad.genomes.r2.0.2.sites.GRCh38.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS,segdup,lcr





vep --input_file Sample_15120+Sample_15121.FINALmutect2.vcf --output_file ./Sample_15120+Sample_15121.FINALmutect2.vep.vcf \
	 --dir /fdb/VEP/101/cache --fasta customDB/GRCh38.d1.vd1.fa \
	 --plugin ExAC,/data/tangw3/VEP/customDB/ExAC.0.3.GRCh38.vcf.gz \
	 --custom /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas,AF_popmax,variant_type,segdup,lcr \
	 --custom /data/tangw3/VEP/customDB/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADe,vcf,exact,0,non_cancer_AF_afr,non_cancer_AF_amr,non_cancer_AF_asj,non_cancer_AF_eas,non_cancer_AF_fin,non_cancer_AF_nfe,non_cancer_AF_oth,non_cancer_AF_sas,non_cancer_AF_popmax,AF_popmax,variant_type,segdup,lcr \
	 --species homo_sapiens --assembly GRCh38 \
	 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 \
	 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf --offline --pubmed \
	 --fork 16 --polyphen b --af --af_1kg --af_esp --af_gnomad --regulatory --force_overwrite --pick

vcf2maf.pl --input-vcf Sample_15120+Sample_15121.FINALmutect2.vep.vcf --output-maf example.maf   \
	--inhibit-vep --ncbi-build GRCh38 --ref-fasta  customDB/GRCh38.d1.vd1.fa \
	--tumor-id Sample_15120 --normal-id Sample_15121  \
	--vep-forks 16 --retain-info CSQ --filter-vcf /data/tangw3/VEP/customDB/ExAC.0.3.GRCh38.vcf.gz 






#####################################################################
######################################################################
###########   make vep_vcf2maf.sh ####################################
#####################################################################




#! /usr/bin/bash

SWARM_FILE="vcf2maf.swarm"
OUT_DIR="maf"
VCF_SUFFIX=".FINALmutect2.vcf"

############# hg38 #############
SPECIES="homo_sapiens"
GENOME="GRCh38" 
GENOME_LOC="/data/tangw3/VEP/customDB/GRCh38.d1.vd1.fa" 
FILTER_SITES="/data/tangw3/VEP/customDB/ExAC.0.3.GRCh38.vcf.gz" 


### check exists ###
SWARM_FILE="vcf2maf.swarm"



if [ -e $SWARM_FILE ]; then
    rm $SWARM_FILE
    touch $SWARM_FILE
fi

if [ ! -e $OUT_DIR ]; then
    mkdir -p $OUT_DIR
fi

### write swarm files ###

for file_name in *.vcf; do    


##### VEP #####

    VEP_out="$OUT_DIR/$(basename -s $VCF_SUFFIX $file_name).vep.vcf"

##### vcf2maf ####
	
      pairs=($(basename -s $VCF_SUFFIX "$file_name" | sed -e s/\+/\\n/))  ## Cleans and splits the filename to extract tumor ID
   tumor_id=${pairs[1]}
    norm_id=${pairs[0]}  

    in_file="$OUT_DIR/$(basename -s $VCF_SUFFIX $file_name).vep.vcf"
   out_file="$OUT_DIR/$(basename -s $VCF_SUFFIX $file_name).maf"

     MY_CMD="date; module load VEP/101; module load vcf2maf/1.6.19; \
     vep --input_file $file_name --output_file $VEP_out \
	 --dir /fdb/VEP/101/cache --fasta $GENOME_LOC \
	 --plugin ExAC,$FILTER_SITES \
	 --custom /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth,AF_sas,AF_popmax,variant_type,segdup,lcr \
	 --custom /data/tangw3/VEP/customDB/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADe,vcf,exact,0,non_cancer_AF_afr,non_cancer_AF_amr,non_cancer_AF_asj,non_cancer_AF_eas,non_cancer_AF_fin,non_cancer_AF_nfe,non_cancer_AF_oth,non_cancer_AF_sas,non_cancer_AF_popmax,AF_popmax,variant_type,segdup,lcr \
	 --species $SPECIES --assembly $GENOME \
	 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 \
	 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --format vcf --offline --pubmed \
	 --fork 16 --polyphen b --af --af_1kg --af_esp --af_gnomad --regulatory --force_overwrite --pick; 
vcf2maf.pl --input-vcf $in_file --output-maf $out_file \
           --inhibit-vep \  ####--vep-path \$VEP_HOME --vep-data \$VEPCACHEDIR if you need VEP automatically
           --species $SPECIES --ncbi-build $GENOME --ref-fasta $GENOME_LOC --filter-vcf $FILTER_SITES \
           --tumor-id "$tumor_id" --normal-id "$norm_id" \  
           --retain-info CSQ \
	   --vep-forks 16; date"
    echo $MY_CMD >> $SWARM_FILE
done

echo "Submitting swarm job..."
swarm_cmd="swarm -f $SWARM_FILE -g 32 --partition quick"
echo $swarm_cmd
#eval $swarm_cmd

echo "Done!"

#####################################################################
######################################################################
###########   make vep_vcf2maf.sh END ####################################
#####################################################################








VCF_FILES=($(ls *$VCF_SUFFIX))

for file_name in "${VCF_FILES[@]}"; do   echo $file_name;done


for file_name in *.vcf;do
      pairs=($(basename -s $VCF_SUFFIX $file_name | sed -e s/\+/\\n/))  ## Cleans and splits the filename to extract tumor ID
   tumor_id=${pairs[0]}
    norm_id=${pairs[1]}  
 echo $tumor_id
 echo $norm_id
done




--appris, --max_af, --var_synonyms,--mane

### VEP config ####

species human
assembly GRCh38
format vcf
force_overwrite
fork 16

cache
dir_cache /data/tangw3/VEP/customDB 
offline
fasta  customDB/GRCh38.d1.vd1.fa
show_cache_info

custom /data/tangw3/VEP/customDB/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AF_NFE,variant_type,segdup

