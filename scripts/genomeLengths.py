#!/usr/bin/env python3

from Bio import SeqIO

import click
import os

"""
Counts total genome length by summing contig lengths
"""

def countGenomeLength(input, prefix:str)->str:
    sum = 0
    for contig in input:
        sum += len(str(contig.seq))
    return f'{prefix}\t{sum}'

@click.command()
@click.option('--input', '-i', help='Specify a genome contig file', type=click.Path(exists=True))

def execute(input):
    prefix = input.split('/')[-1].split('.')[0]
    with open(input, 'r') as handle:
        handle = SeqIO.parse(handle, 'fasta')
        print(type(handle))
        result = countGenomeLength(handle, prefix)


if __name__ == '__main__':
    execute()
