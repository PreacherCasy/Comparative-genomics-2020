#!/usr/bin/env python3

import click
import os
from Bio import SeqIO
from Bio.Seq import Seq

"""
Translate all the nucleotide sequences from the file, then write the protein entries to the output file
"""


def translate(input:str, output:str):

    output_list = list()
    with open(input, 'r') as handle:
        for entry in SeqIO.parse(handle, 'fasta'):
            entry.seq = entry.seq.translate()
            output_list.append(entry)

    with open(output, 'w') as handle:
        SeqIO.write(output_list, handle, 'fasta')


@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path(exists=False))

def execute(input, output):
    translate(input, output)


if __name__ == '__main__':
    execute()
