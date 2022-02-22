import argparse
import glob
import pandas as pd

BLAST_HEADER ='query target match/tLen mismatch gapOpen qstart qend tstart tend score eval'.split()

parser = argparse.ArgumentParser(description='Script to convert HHsuite results to tsv')

# Define parser arguments
parser.add_argument('result_folder', type=str, help='Folder containing the HHblits search results')
parser.add_argument('output_path', type=str, help='Path to write parsed report to')
parser.add_argument('phrog_file', type=str, help='Path to Phrog annotation table')
parser.add_argument('prodigal_file', type=str, help='Path to prodigal coordinate table')

# Get input arguments
args = parser.parse_args()
result_folder = args.result_folder
output_path = args.output_path
phrog_file = args.phrog_file
prod_file = args.prodigal_file

# Load annotation table
phrog_anno = pd.read_csv(phrog_file, delimiter = '\t')

# Parse prodigal coordinates
prodigal_table = []
with open(prod_file,'r') as prodigal_file:
    for line in prodigal_file:
        if 'seqhdr=' in line:
            contig_id = line.split('seqhdr="')[-1].rstrip('"\n')

        if not line.startswith('#'):
            orf_id, start, end, strand = line.strip().lstrip('>').split('_')
            prodigal_table.append({
                'query' : contig_id + '_' + orf_id,
                'start' : int(start),
                'end': int(end),
                'strand': strand
            })
prodigal_table = pd.DataFrame(prodigal_table)


header = ['protein_id'] + BLAST_HEADER
hmm_table = []
total, no_match = 0, 0
for file in glob.glob(result_folder.rstrip('/') + '/*.txt'):
    total +=1
    with open(file) as in_file:
        result = in_file.readline().strip()
        if result:
            data = dict(zip(BLAST_HEADER, result.split('\t')))
            data['phrog'] = int(data['target'].split('_')[-1])
            hmm_table.append(data)
        else:
            no_match +=1
print('[INFO] {} out of the {} proteins did not have a match'.format(no_match, total))

# Convert to pandas dataframe and left join annotation table
hmm_table = pd.DataFrame(hmm_table)
final_table = hmm_table.merge(phrog_anno, on = ['phrog'], how='left')
final_table = final_table.merge(prodigal_table, on = ['query'], how = 'left')
print(final_table.head())

# Write to output file
final_table.to_csv(output_path, '\t', index = False)

