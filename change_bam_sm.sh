samtools view -H $BAM | sed "s/SM:Oldname/SM:Newname/" | samtools reheader - $BAM
