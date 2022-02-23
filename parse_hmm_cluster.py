import argparse
import glob
import pandas as pd


BLAST_HEADER ='query target match/tLen mismatch gapOpen qstart qend tstart tend score eval'.split()

parser = argparse.ArgumentParser(description='Script to convert HHsuite results to tsv')
#
# Define parser arguments
parser.add_argument('result_folder', type=str, help='Folder containing the HHblits search results')
parser.add_argument('output_path', type=str, help='Path to write parsed report to')
parser.add_argument('phrog_file', type=str, help='Path to Phrog annotation table')
parser.add_argument('prodigal_file', type=str, help='Path to prodigal coordinate table')
parser.add_argument('cluster_tsv', type=str, help='Path to cluster.tsv file')
parser.add_argument('cluster_stats', type=str, help='Path to cluster stat file')

# Get input arguments
args = parser.parse_args()
result_folder = args.result_folder
output_path = args.output_path
phrog_file = args.phrog_file
prod_file = args.prodigal_file
cluster_file = args.cluster_tsv
cluster_stat_file = args.cluster_stats

# Load annotation table
phrog_anno = pd.read_csv(phrog_file, delimiter = '\t')

# Parse prodigal coordinates
# We can join this with the
prodigal_table = []
with open(prod_file,'r') as prodigal_file:
    for line in prodigal_file:
        if 'seqhdr=' in line:
            contig_id = line.split('seqhdr="')[-1].rstrip('"\n')

        if not line.startswith('#'):
            orf_id, start, end, strand = line.strip().lstrip('>').split('_')
            prodigal_table.append({
                'member' : contig_id + '_' + orf_id,
                'start' : int(start),
                'end': int(end),
                'strand': strand
            })
prodigal_table = pd.DataFrame(prodigal_table)

# Read the cluster file (no header)
clusters = pd.read_csv(cluster_file, header = None, delimiter = '\t')
clusters.columns = ['reference', 'member']

# Read the cluster stat file
cluster_stat = pd.read_csv(cluster_stat_file, delimiter = '\t')
pd.set_option('display.expand_frame_repr', False)

# Join the cluster tables
clusters = pd.merge(clusters, cluster_stat, on = 'reference', how = 'left')

# Parse the single result files
header = ['cluster_id'] + BLAST_HEADER
rows = []
out_base = result_folder.rstrip('/')
for file in glob.glob(out_base + '/*.txt'):
    cluster_id = file.split('/')[-1].rstrip('.txt')
    with open(file) as in_file:
        data = dict(zip(BLAST_HEADER, in_file.readline().split('\t')))
        data['cluster_id'] = cluster_id
        print(cluster_id, data)
        rows.append(data)
search_results = pd.DataFrame(rows)

# Combine the clusters and search results
search_results = search_results.astype({'cluster_id': 'int64'})
final_table = pd.merge(search_results, clusters, on = 'cluster_id', how = 'left')

# # Also add the annotation
phrog_anno['target'] = ['phrog_{}'.format(x) for x in phrog_anno['phrog']]
final_table = pd.merge(final_table, phrog_anno, how = 'left')
final_table.to_csv('output_test.txt', sep = '\t', index = False)

# Also add the gene coordinate annotation
final_table = pd.merge(final_table, prodigal_table, how = 'left')

# Save to output path
final_table.to_csv(output_path, index = False, sep = '\t')


