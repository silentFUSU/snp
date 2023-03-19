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
' | samtools view -b -o S04_primary-atac_new.bam -
