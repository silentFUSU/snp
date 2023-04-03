java -jar ~/software/picard.jar AddOrReplaceReadGroups I=${i}.realn.bam O=${i}.bam RGID=${i} RGLB=${i} RGPL=illumina RGSM=${i} RGPU=unit1
