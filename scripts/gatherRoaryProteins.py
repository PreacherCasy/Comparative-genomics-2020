#!/usr/bin/env python3

from Bio import SeqIO
from itertools import combinations

import click
import os
import pandas as pd
import subprocess

"""
Collects the core genes based on the Roary gene presence/absence matrix
"""

def extractUniquesRoary(input:str):

    """
    Extracts the names of clusters containing exactly one ortholog per genome
    """

    prmat = pd.read_csv(input, sep=',', header=0, index_col=False)
    gencols = [x for x in prmat.columns.values if x.find('GCF') > -1]
    genes = [prmat.iloc[x,0] for x in range(0, prmat.shape[0]) if int(prmat.iloc[x,3]) == len(gencols) and int(prmat.iloc[x,4]) == len(gencols)]

    return genes

def checkCollinearity(input:list)->bool:
    """"
    A very straightforward function to verify that the genes from the same genome share the same index across the clusters
    """
    master_list = []
    for file in input:
        cluster_list = []
        with open(file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                cluster_list.append(record.id.split('_')[0])
        master_list.append(cluster_list)

    for a, b in combinations(master_list, 2):
        if a != b:
            return False

    return True


@click.command()
#@click.option('--matrix', '-m', help='A Roary-derived distance matrix', type=click.Path(exists=True))
@click.option('--input', '-i', help='A directory containing Roary output', type=click.Path(exists=True))
@click.option('--output', '-o', help='A directory to store the gathered sequences', type=click.Path(exists=False))
@click.option('--sanity', '-S', help='Perform file collinearity check', is_flag=True)

def execute(input:str, output:str, sanity:bool):
    if not os.path.exists(output):
        os.makedirs(output)

    matr = os.path.join(input, 'gene_presence_absence.csv')
    genes = extractUniquesRoary(matr)

    refs = os.path.join(input, 'pan_genome_sequences')
    for file in os.listdir(refs):
        prefix = file.split('.')[0]
        if prefix in genes:
            init_pos = os.path.join(refs, file)
            fin_pos = os.path.join(output, file)
            subprocess.run(f'cp {init_pos} {fin_pos}', shell=True)

    if sanity:
        files = list(map(lambda x: os.path.join(output, x), os.listdir(output)))
        collinearity = checkCollinearity(files)
        if collinearity:
            print('The sequence order is preserved')
        else:
            raise Exception('Sequence order is not preserved throughout the clusters!')

if __name__ == '__main__':
    execute()
