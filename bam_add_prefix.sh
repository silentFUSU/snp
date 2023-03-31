#add prefix to barcode

func() {
    echo "Usage:"
    echo "add_prefix [-f bam_dir] [-p prefix] [-o out_dir]"
}

while getopts ':h:f:p:o:' OPT; do
    case $OPT in
        f) bam_dir=${OPTARG};;
        p) prefix="$OPTARG";;
        o) out_dir=${OPTARG};;
        h) func;;
        ?) func;;
    esac
done

echo "f:"$bam_dir
echo "p:"${prefix}_
echo "o:"$out_dir
# samtools view -h ${bam_dir}  
samtools view -h ${bam_dir} | awk -v OFS="\t" -v prefix=${prefix}_ '
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
' | samtools view -b -o ${out_dir} -
samtools index ${out_dir}
