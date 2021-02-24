

#! /usr/bin/bash
#wget -c https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/copyCat/GRCh38/gaps.bed
#cat gaps.bed | awk '{print "chr"$0}' > tmp && mv tmp gaps.bed
#cnvkit.py access /data/tangw3/reference/hg38.fa -x gaps.bed -o access-excludes.hg38.bed

#cnvkit.py batch /data/WES_breastcancer/analysis/pipeliner_output/wes_hg38/Sample_2.recal.bam --normal /data/WES_breastcancer/analysis/pipeliner_output/wes_hg38/Sample_1.recal.bam \
#	--targets /data/tangw3/reference/cnvkit/S07604624_Padded.bed --annotate /data/tangw3/reference/cnvkit/refFlat.txt --fasta /data/tangw3/reference/hg38.fa --access /data/tangw3/reference/cnvkit/access-excludes.hg38.bed \
#	--output-reference /data/tangw3/reference/cnvkit/test/my_reference.cnn --output-dir /data/tangw3/reference/cnvkit/test/ \
#	--method hybrid -p 16 --diagram --scatter --drop-low-coverage


###################################
##### input and output params #####
###################################

normal_DIR="/data/tangw3/cnvkit/normal"
tumor_DIR="/data/tangw3/cnvkit/tumor"

OUT_DIR="/data/tangw3/cnvkit/results"

target_file="/data/tangw3/reference/cnvkit/S07604624_Padded.bed"
annotate_file="/data/tangw3/reference/cnvkit/refFlat.txt"
fasta_file="/data/tangw3/reference/hg38.fa"
access_file="/data/tangw3/reference/cnvkit/access-excludes.hg38.bed"


module load cnvkit

if [ ! -e $OUT_DIR ]; then
    mkdir -p $OUT_DIR
fi

cnvkit.py batch $tumor_DIR/*.bam --normal $normal_DIR/*.bam \
    --targets $target_file --annotate $annotate_file \
    --fasta $fasta_file --access $access_file \
    --output-reference $OUT_DIR/reference.cnn --output-dir $OUT_DIR \
    --method hybrid -p 32 --diagram --scatter --drop-low-coverage






############## gistic2
#!/bin/sh
set -e

module load gistic

basedir=gistic2_results_0.3
mkdir gistic2_results_0.3

segfile=gistic.segments
markersfile=$GISTIC_HOME/examplefiles/markersfile.txt
refgenefile=$GISTIC_HOME/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
alf=$GISTIC_HOME/examplefiles/arraylistfile.txt
cnvfile=CNV.hg38.bypos.111213.txt

## call script that sets MCR environment and calls GISTIC executable 
gistic2 -b $basedir -seg $segfile -refgene $refgenefile -cnv $cnvfile -genegistic 1 -qvt 0.05 -js 20 -ta 0.3 -td 0.3 -smallmem 0 -broad 1 -brlen 0.98 -conf 0.99 -res 0.01 -armpeel 1 -savegene 1 -gcm extreme

mv gistic2_results_0.3 gistic2_results_0.3.js20.res.0.01.conf0.99








