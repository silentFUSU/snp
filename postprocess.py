import pysam

bamfile = pysam.AlignmentFile('bam.bam', 'rb')
cell_counts = {}
for read in bamfile.fetch('chr21', 8214763, 8214764):
    if read.query_sequence[read.query_alignment_start] == 'G':
        if read.has_tag('CB'):
            cell_id = read.get_tag('CB')
            cell_counts[cell_id] = cell_counts.get(cell_id, 0) + 1

bamfile.close()

with open('cell_ids_rec.txt', 'w') as f:
    for cell_id, count in cell_counts.items():
        f.write(cell_id + '\n')
