samtools view -h S04_primary-atac.sam | awk -v OFS="\t" -v prefix="pri_" '
    BEGIN { FS = "\t" } 
    $1 ~ /^@/ { print $0 } 
    $1 !~ /^@/ { 
        for (i=1; i<=NF; i++) {
            if ($i ~ /^CB:Z:/) {
                sub(/^CB:Z:/, "CB:Z:"prefix, $i);
            }
        }
        print $0;
    }
' | samtools view -b -o S04-primary-atac_new.bam -
samtools index S04-primary-atac_new.bam

samtools view -h S04-recurrent-atac.bam | awk -v OFS="\t" -v prefix="rec_" '
    BEGIN { FS = "\t" } 
    $1 ~ /^@/ { print $0 } 
    $1 !~ /^@/ { 
        for (i=1; i<=NF; i++) {
            if ($i ~ /^CB:Z:/) {
                sub(/^CB:Z:/, "CB:Z:"prefix, $i);
            }
        }
        print $0;
    }
' | samtools view -b -o S04-recurrent-atac_new.bam -
samtools index S04-recurrent-atac_new.bam

samtools merge S04.bam S04-primary-atac_new.bam S04-recurrent-atac_new.bam
samtools index S04.bam
