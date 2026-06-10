#!/bin/bash

for sample in p5_s1 p5_s2 p5_s3 p15_s1 p15_s2 p15_s3; do
    (
    start=$(date +%s)
    STAR --runThreadN 1 --genomeDir ./ref/star/ --readFilesIn \
            ./rawdata/1.Fastq/hADSC_${sample}_r1.fastq.gz \
            ./rawdata/1.Fastq/hADSC_${sample}_r2.fastq.gz \
        --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ./res/star/hADSC_${sample}_ --runMode alignReads \
        > ./log/align_star_${sample}.log 2>&1 \
        && samtools index ./res/star/hADSC_${sample}_Aligned.sortedByCoord.out.bam
    end=$(date +%s)
    echo "[$(date)] ${sample} done | elapsed: $((end - start))s" >> ./log/align_star_${sample}.log
    ) &
done
wait
echo "[$(date)] STAR complete"
