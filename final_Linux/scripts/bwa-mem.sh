#!/bin/bash

for sample in p5_s1 p5_s2 p5_s3 p15_s1 p15_s2 p15_s3; do
    (
    start=$(date +%s)
    bwa mem -t 1 ./ref/bwa/GRCh38p14 ./rawdata/1.Fastq/hADSC_${sample}_r1.fastq.gz ./rawdata/1.Fastq/hADSC_${sample}_r2.fastq.gz \
	2> ./log/align_bwa_${sample}.log \
        | samtools sort -@ 1 -T ./res/bwa/${sample}_tmp -o ./res/bwa/hADSC_${sample}_sorted.bam \
        && samtools index ./res/bwa/hADSC_${sample}_sorted.bam
    end=$(date +%s)
    echo "[$(date)] ${sample} done | elapsed: $((end - start))s" >> ./log/align_bwa_${sample}.log
    ) &
done
wait
echo "[$(date)] BWA-MEM complete"
