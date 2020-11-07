#! /usr/bin/env python3

#from alignClusters import correctDuplicateNames
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo import TreeConstruction

import click
import numpy as np
import pandas as pd
""""
Calculates distance matrix for an alignment and writes it to the TSV file
"""

def openAlignment(input:str, format:str)->MultipleSeqAlignment:
    with open(input, 'r') as handle:
        alignment = AlignIO.read(handle, format)

        return alignment


def distMat(alignment:MultipleSeqAlignment, method:str):
    calculator = TreeConstruction.DistanceCalculator('identity')
    distance = calculator.get_distance(alignment)

    return distance

@click.command()
@click.option('--input', '-i', help='Specify an input alignment file', type=str)
@click.option('--output', '-o', help='Specify an output table file', type=str)
@click.option('--format', '-f', help='Specify an alignment format (default is "fasta"', type=str, default='fasta')
@click.option('--method', '-m', help='Set a distance evaluation method (default is "identity")', type=str, default='identity')

def execute(input, output, format, method):
    """
    Opens an alignment file, calculates a distance matrix, then writes it as a 'tsv' dataframe
    """

    alignment = openAlignment(input, format)
    distmat = distMat(alignment, method)
    matrix = np.zeros([len(distmat), len(distmat)])
    for i in range(0, len(matrix)):
        for j in range(0, i+1):
            matrix[i][j] = 1 - distmat[i][j]
    matrix = matrix + matrix.T
    df = pd.DataFrame(matrix, columns = distmat.names, index = distmat.names)
    df.to_csv(output, header=True, index=True, sep='\t')

if __name__ == '__main__':
    execute()
