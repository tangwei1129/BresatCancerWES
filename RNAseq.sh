#! /bin/bash

module load cutadapt
module load samtools
module load STAR
module load rsem
module load picard


adapter=/data/tangw3/reference/adapter.fa
GTF=/data/tangw3/reference/rsem_ref/Rattus_norvegicus.Rnor_6.0.101.chr.gtf
refFlat=/data/tangw3/reference/rsem_ref/refFlat.txt
star_genome=/fdb/STAR_current/Ensembl/Rnor_6.0/genes-50/
rsem_ref=/data/tangw3/reference/rsem_ref/rat_ref 
numthread=8
samplename=sample

cat *_R1_001.fastq.gz *_R1.fastq.gz > input_R1.fastq.gz
cat *_R2_001.fastq.gz *_R2.fastq.gz > input_R2.fastq.gz

cutadapt \
    -j $numthread \
    -b file:$adapter -B file:$adapter \
    --trim-n -n 5 -O 5 -q 10,10 -m 35:35 \
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz input_R1.fastq.gz input_R2.fastq.gz

STAR \
    --runThreadN $numthread \
    --genomeDir $star_genome \
    --readFilesIn trimmed_R1.fastq.gz trimmed_R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM \
    --sjdbGTFfile $GTF

rsem-calculate-expression -bam --paired-end --estimate-rspd Aligned.toTranscriptome.out.bam $rsem_ref $samplename

java -Xmx32g -XX:ParallelGCThreads=8 -jar $PICARDJARPATH/picard.jar CollectRnaSeqMetrics \
    REF_FLAT=$refFlat \
    INPUT=Aligned.out.bam \
    OUTPUT= RnaSeqMetrics.txt \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND VALIDATION_STRINGENCY=LENIENT

samtools sort Aligned.out.bam > Aligned.out.sorted.bam

java -Xmx32g -XX:ParallelGCThreads=8 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
    INPUT=Aligned.out.sorted.bam \
    OUTPUT= MKDUP.bam \
    METRICS_FILE=sample.bam.metric \
    ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT

rm RNAseq.sh
