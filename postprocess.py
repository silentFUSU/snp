import pysam

bamfile = pysam.AlignmentFile('bam.bam', 'rb')

with open("snp_position.txt", "r") as f:
    for line in f:
        fields = line.strip().split()
        mutation_name = "-".join(fields[:4])
        cell_counts = {}
        for read in bamfile.fetch(fields[0],int(fields[1])-1,int(fields[1])):
            if read.query_sequence[read.query_alignment_start] == fields[2]:
                if read.has_tag('CB'):
                    cell_id = read.get_tag('CB')
                    cell_counts[cell_id] = cell_counts.get(cell_id, 0) + 1
        with open(f'{mutation_name}-pri1.txt', 'w') as f:
            for cell_id, count in cell_counts.items():
                f.write(cell_id + '\n') 
        cell_counts = {}
        for read in bamfile.fetch(fields[0],int(fields[1])-1,int(fields[1])):
            if read.query_sequence[read.query_alignment_start] == fields[3]:
                if read.has_tag('CB'):
                    cell_id = read.get_tag('CB')
                    cell_counts[cell_id] = cell_counts.get(cell_id, 0) + 1
        with open(f'{mutation_name}-pri2.txt', 'w') as f:
            for cell_id, count in cell_counts.items():
                f.write(cell_id + '\n') 
bamfile.close()
