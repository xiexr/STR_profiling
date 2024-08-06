echo ========step 1: run alignment by bwa========
bwa mem -M -t 20 -R "@RG\\tID:$samplename\\tLB:$samplename\\tSM:$samplename\\tPU:$samplename\\tPL:Illumina" /home/dragon/DB/riceIRGSP/genome/riceIRGSP.fa $sample_R1.fq.gz $sample_R2.fq.gz | samtools view -F 4 -Sb - | samtools sort -@ 20 -O bam -o $sample.sort.bam
echo ========step 2: remove duplication ==========
samtools rmdup $samplename.srt.bam $samplename.rp.bam 
samtools index $samplename.rp.bam