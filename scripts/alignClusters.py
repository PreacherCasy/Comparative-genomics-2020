#!/usr/bin/env python3

from Bio import Align
from Bio import AlignIO
from Bio import SeqIO
from Bio.Phylo import TreeConstruction

import click
import os
import subprocess
import numpy as np
import pandas as pd

"""
Aligns sequences with aligner of your choice, then produces a distance matrix for each alignent and a pivot table with minimum and maximum identity value for each distance matrix
"""

def readNames(input:str)->list:
    """
    Read a single-column list of single-copy orthologs
    """

    output = []
    with open(input, 'r') as handle:
        for line in handle.readlines():
            output.append(line.rstrip('\n'))

    return output


def alignFilewise(input:str, outdir:str, aligner:str, *args):
    """
    Align sequences using either MAFFT or PRANK. The default option is MAFFT. Note that the aligner should be reachable from the PATH
    """

    outdir = os.path.abspath(outdir)
    output = f'{input.split("/")[-1].split(".")[0]}.aln'

    if aligner == 'mafft':
        command = [aligner] + list(args) + [f'{input}', '>', f'{os.path.join(outdir, output)}']
    elif aligner == 'prank':
        command = [aligner] + list(args) + [f'-d={input}', f'-o={os.path.join(outdir, output)}']
    subprocess.run(' '.join(command), shell=True, check=True)


def correctDuplicateNames(input:Align.MultipleSeqAlignment):
    """
    Corrects duplicate sequence names in multisequence alignments. Might be useful when dealing with independent genome annotation outputs
    """

    ids = []
    for i in range(0, len(input._records)):
        prefix = input._records[i].id.split('.')[0]
        ids.append(prefix)
        if ids.count(prefix) > 1:
            input._records[i].id = f'{prefix}.{ids.count(prefix) + 1}'
            input._records[i].name = input._records[i].id

    return input

def calculateDistanceMatrix(input:str)->pd.DataFrame:
    """
    Calculates distance matrix for a given alignment file in fasta format
    """

    # open an alignment file
    with open(input, 'r') as handle:
        aln = correctDuplicateNames(AlignIO.read(handle, 'fasta'))

    # calculate a distance matrix
    calc = TreeConstruction.DistanceCalculator(model='identity')
    distance = calc.get_distance(aln)

    # generate a symmetrical distance matrix
    matr = np.zeros([len(distance.matrix), len(distance.matrix)])
    for i in range(0, len(matr)):
        for j in range(0, i+1):
            matr[i][j] = 1 - distance.matrix[i][j]
    matr = matr + matr.T - np.diag(matr.diagonal())


    # generate a data frame
    df = pd.DataFrame(matr, columns = distance.names, index = distance.names)

    return df


@click.command()
@click.option('--names', '-n', help='A file containing a list of cluster names to use', type=click.Path(exists=True)) # a name file
@click.option('--input', '-i', help='An input directory containing orthologous sequences', type=click.Path(exists=True))
@click.option('--output', '-o', help='An output directory', type=click.Path(exists=False))
# @click.option('--pivot', '-p') # a path to a pivot table containing minimum and maximum distance values per each sequence
@click.option('--aligner', '-a', help='Specify whether PRANK or MAFFT should be used for MSA. Note that it should be reachable from the current $PATH', type=click.Choice(['prank', 'mafft'], case_sensitive=False), default='mafft')
@click.option('--opts', '-m', help='Additional parameters for the aligner chosen',  multiple=True, type=str, default='')


def execute(names, input, output, aligner, opts):
    """
    Aligns all sequences found in the input directory filewise, then returns distance matrices for each alignment and a pivot table containing minimum and maximum identity values for each orthogroup
    """

    if not os.path.exists(output):
        os.makedirs(output)

    name_list = readNames(names)
    pivot_dict = dict()
    files = os.listdir(input)
    for file in files:
        prefix = file.split('.')[0]
        if prefix in name_list:
            if not os.path.exists(os.path.join(output, 'alignments')):
                os.makedirs(os.path.join(output, 'alignments'))
            print(os.path.join(output, 'alignments'))
            alignFilewise(os.path.join(input, file), os.path.join(output, 'alignments'), aligner, *mafft_opts)
            dm = calculateDistanceMatrix(os.path.join(output, 'alignments', f'{prefix}.aln'))
            if not os.path.exists(os.path.join(output, 'distmats')):
                os.makedirs(os.path.join(output, 'distmats'))
            with open(os.path.join(output, 'distmats', f'{prefix}_distmat.tsv'), 'w') as handle:
                dm.to_csv(handle, header=True, index=True, sep='\t')
            matr = dm.to_numpy()
            matr[matr == 1] = 0
            pivot_dict[prefix] = [np.min(matr[np.nonzero(matr)]),
                                  np.max(matr[np.nonzero(matr)])]

    pivot_df = pd.DataFrame({'Orthogroup': list(pivot_dict.keys()),
                             'MinIdent': list(map(lambda x: x[0], pivot_dict.values())),
                             'MaxIdent': list(map(lambda x: x[1], pivot_dict.values()))})
    with open(os.path.join(output, 'summary.tsv'), 'w') as handle:
        pivot_df.to_csv(handle, header=True, index=False, sep='\t')


if __name__ == '__main__':
    execute()

