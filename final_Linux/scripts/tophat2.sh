#!/bin/bash

for sample in p5_s1 p5_s2 p5_s3 p15_s1 p15_s2 p15_s3; do
    (
    start=$(date +%s)
    tophat2 -p 1 -G ./ref/gencode.v49.annotation.gtf -o ./res/tophat2/${sample} ./ref/tophat2/GRCh38p14 \
        ./rawdata/1.Fastq/hADSC_${sample}_r1.fastq.gz \
        ./rawdata/1.Fastq/hADSC_${sample}_r2.fastq.gz \
        2> ./log/align_tophat2_${sample}.log
    end=$(date +%s)
    echo "[$(date)] ${sample} done | elapsed: $((end - start))s" >> ./log/align_tophat2_${sample}.log
    ) &
done
wait
echo "[$(date)] TopHat2 complete"
