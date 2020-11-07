#! /usr/bin/env python3

import click
from Bio import SeqIO

def checkLengths(file):

    with open(file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            print(f'{file.split("/")[-1]}\t{record.id}\t{len(record.seq)}')


@click.command()
@click.argument('input', type=click.Path(exists=True))

def execute(input):
    checkLengths(input)


if __name__ == '__main__':
    execute()

