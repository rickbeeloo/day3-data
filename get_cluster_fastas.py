from pyfaidx import Fasta
import os
import sys
import argparse
from tabulate import tabulate

parser = argparse.ArgumentParser(description='Script to generate a FASTA file for each cluster '
                                             'containing at least 1 prodigal protein and 1 RefSeq protein')


# Command line arguments
parser.add_argument('output_folder', type=str, help='Path to the output folder')
parser.add_argument('proteins', type=str, help='Path to prodigal protein FASTA (.faa) file')
parser.add_argument('clusters', type=str, help='Path to mmseqs2 cluster tsv')
parser.add_argument('--min_seqs', type=int,  help='Minimum number of sequences that should be present'
                                                  ' in a cluster to include it (DEFAULT = 2)', default=argparse.SUPPRESS)

# Parse the provided command line arguments
args = parser.parse_args()
out_folder = args.output_folder
prodigal_proteins = args.proteins  #'prodigal/l200.proteins_with_refseq.faa'
cluster_file = args.clusters # 'clusters_with_refseq.tsv'
if hasattr(args, 'min_seqs'):
    min_seqs = args.min_seqs
else:
    print('Min seqs not provided')
    min_seqs = 2

# Provide the user with an overview of the input
used_arguments = [['output folder', out_folder], ['prodigal proteins', prodigal_proteins], ['cluster file', cluster_file], ['min sequences', min_seqs]]
print(tabulate(used_arguments, headers=["Argument", "Value"]) + '\n')

# Create output directory
# - needs to be empty to work with HH-suite so
# - we abort when it's not
def create_out_dir(path):
    try:
        os.mkdir(path)
    except OSError as error:
        # Directory exists, so check if it's empty
        if not len(os.listdir(path)) == 0:
            print('[ERROR] Directory is not empty!')
            sys.exit(0)


# Construct output paths
fasta_out = out_folder.rstrip('/') + '/Fasta'
msa_out = out_folder.rstrip('/') + '/MSA'
stat_out = out_folder.rstrip('/') + '/Scores'

# Create output folders
create_out_dir(out_folder)
create_out_dir(fasta_out)
create_out_dir(msa_out)
create_out_dir(stat_out)

# Function to write a cluster to an output file
def write_cluster(cluster_id, fasta_out_path, fasta):
    with open(fasta_out_path, 'w') as out:
        out.write(fasta)

# Function to check if a cluster contains a by prodigal predicted protein
def prodigal_count(accs):
    return sum(['NODE' in acc for acc in accs])

# Each cluster is associated with an "unique" reference protein id, but it's a bit
# messy to use these as out MSA outputs. So we map the protein ids to unique numbers
cluster_id = 0
cluster_map = dict()

# We create a FASTA index to be able to quickly grab sequences from
# the prodigal protein file
prodigal_proteins = Fasta(prodigal_proteins)

# Now we parse the cluster file
last_ref = None
cluster_fasta = ''
member_ids = []

wrote = 0

with open(cluster_file,'r') as in_file, open(out_folder.rstrip('/') + '/cluster_stats.txt','w') as out:
    # Write header
    out.write('cluster_id\treference\tmembers\tprodigal\n')

    for line in in_file:

        # Get the (cluster) reference protein id and the member protein id
        ref, member = line.rstrip().split('\t')

        # Get the sequence from the index and create a FASTA line
        seq = prodigal_proteins[member]
        fasta = '>{}\n{}\n'.format(member, seq)

        # Check if we should append to the current cluster or create a
        # new cluster
        if ref != last_ref:

            # Write current cluster if it has at least one prodigal identifier
            prodigal_ones = prodigal_count(member_ids)
            if prodigal_ones >= 1 and len(member_ids) >= min_seqs:
                wrote +=1
                fasta_path = fasta_out + '/{}.fasta'.format(cluster_id)
                out.write('\t'.join(map(str, [cluster_id, cluster_map[cluster_id], len(member_ids), prodigal_ones])) + '\n')
                write_cluster(cluster_id, fasta_path, cluster_fasta)


            # Reset vars
            cluster_id += 1
            cluster_map[cluster_id] = ref
            cluster_fasta = ''
            last_ref = ref
            member_ids = []

        cluster_fasta += fasta
        member_ids.append(member)

    # We finished the loop so there might be one cluster left, write that too
    prodigal_ones = prodigal_count(member_ids)
    if prodigal_ones >= 1 and len(member_ids) > min_seqs:
        print(cluster_id, cluster_map[cluster_id])
        out.write('\t'.join(map(str, [cluster_id, cluster_map[cluster_id], len(member_ids), prodigal_ones])) + '\n')
        write_cluster(cluster_id, fasta_path, cluster_fasta)

    print('Done! wrote {} Fasta files'.format(wrote))





