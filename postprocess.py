import pysam
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument("-bf","--bamfile" , help="the path to a bam file")
parser.add_argument("-tf","--txtfile" , help="the path to the snp txt file")
parser.add_argument("-sample", help="sample name")
parser.add_argument("-condition", help="sample condition")

args = parser.parse_args()
bamfile = pysam.AlignmentFile('{}'.format(args.bamfile), 'rb')
result_path='/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/result/merge_data/before_doubletfinder/snp'
outputdir='{}/snp_{}_{}'.format(result_path,args.sample,args.condition)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

#1 ref 2 alt
with open('{}'.format(args.txtfile), "r") as f:
    for line in f:
        fields = line.strip().split()
        mutation_name = "-".join(fields[:4])
        cell_counts = {}
        for read in bamfile.fetch(fields[0],int(fields[1])-1,int(fields[1])):
            if read.query_sequence[read.query_alignment_start] == fields[2]:
                if read.has_tag('CB'):
                    cell_id = read.get_tag('CB')
                    cell_counts[cell_id] = cell_counts.get(cell_id, 0) + 1
        with open('{}/snp_{}_{}/{}-{}1.txt'.format(result_path,args.sample,args.condition,mutation_name,args.condition), 'w') as f:
            for cell_id, count in cell_counts.items():
                f.write(cell_id + '\n') 
        cell_counts = {}
        for read in bamfile.fetch(fields[0],int(fields[1])-1,int(fields[1])):
            if read.query_sequence[read.query_alignment_start] == fields[3]:
                if read.has_tag('CB'):
                    cell_id = read.get_tag('CB')
                    cell_counts[cell_id] = cell_counts.get(cell_id, 0) + 1
        with open('{}/snp_{}_{}/{}-{}2.txt'.format(result_path,args.sample,args.condition,mutation_name,args.condition), 'w') as f:
            for cell_id, count in cell_counts.items():
                f.write(cell_id + '\n') 

bamfile.close()
