#!/bin/bash

for sample in p5_s1 p5_s2 p5_s3 p15_s1 p15_s2 p15_s3; do
    (
    start=$(date +%s)
    hisat2 -p 1 --dta -x ./ref/hisat2/GRCh38p14 \
        -1 ./rawdata/1.Fastq/hADSC_${sample}_r1.fastq.gz \
        -2 ./rawdata/1.Fastq/hADSC_${sample}_r2.fastq.gz \
        2> ./log/align_hisat2_${sample}.log \
        | samtools sort -@ 1 -T ./res/hisat2/${sample}_tmp -o ./res/hisat2/hADSC_${sample}_sorted.bam \
        && samtools index ./res/hisat2/hADSC_${sample}_sorted.bam
    end=$(date +%s)    
    echo "[$(date)] ${sample} done | elapsed: $((end - start))s" >> ./log/align_hisat2_${sample}.log
	) &
done
wait
echo "[$(date)] HISAT2 complete"
