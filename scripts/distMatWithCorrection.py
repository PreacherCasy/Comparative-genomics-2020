#!/usr/bin/env python3

from distMat import distMat, openAlignment
from Bio import Align

import click
import numpy as np
import os
import pandas as pd

"""
Corrects RNA to DNA. An auxiliary script to test whether distance calculator treats DNA and RNA alphabets similarly when calculating identity-based distances.
"""

def correctToDNA(alignment):
    from Bio.Seq import Seq
    for i in alignment:
        i.seq = Seq(str(i.seq).replace('U', 'T'))
    return alignment

@click.command()
@click.option('--input', '-i', help='Specify an input alignment file', type=click.Path(exists=True))
@click.option('--output', '-o', help='Specify an output file', type=str)
@click.option('--format', '-f', help='Specify an alignment format (default is "fasta")', type=str, default='fasta')
@click.option('--method', '-m', help='Set a distance evaluation method (default is "identity")', type=str, default='identity')

def execute(input, output, format, method):
    alignment = openAlignment(input, format)
    alignment = correctToDNA(alignment)
    distance = distMat(alignment, method)
    matrix = np.zeros([distance.__len__(), distance.__len__()])
    for i in range(0, len(matrix)):
        for j in range(0, i + 1):
           matrix[i][j] = 1 - distance[i][j]
    matrix = matrix + matrix.T - np.diag(matrix.diagonal())
    df = pd.DataFrame(matrix, columns = distance.names, index=distance.names)
    df.to_csv(output, header=True, index=True, sep='\t')

if __name__ == '__main__':
    execute()
