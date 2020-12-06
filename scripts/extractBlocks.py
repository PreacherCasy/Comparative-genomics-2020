#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

import click
import _io
import os
import pandas as pd

def applyToAllFastas():
    pass


def parseSibeliazBlocks(input:_io.TextIOWrapper):
    pass


def parseParebrickBlocks(input:_io.TextIOWrapper)->defaultdict:
    output_dict = defaultdict(lambda: defaultdict(list))
    while input:
        try:
            line = next(input).rstrip('\n')
            if line.find('>') > -1:
                block = line.lstrip('>').rstrip('\n')
            else:
                if len(line) > 0:
                    genome, coords = line.split(':')
                    coords, strand = coords.split(' ')
                    start, stop = coords.split('-')
                    output_dict[genome][block] = [start, stop, strand]
        except StopIteration:
            break
    return output_dict


def extractBlocks(file:str, block_dict:defaultdict)->list:

    """
    Extracts blocks contained in the 'block_dict' dictionary.
    Note that the prefix used will contain all the symbols in the filename going before the extention.
    If the resulting prefices do not match the dictionary keys, consider fixing the filenames and/or input block lists manually.
    """

    prefix = '.'.join(file.split('/')[-1].split('.')[0:-1])
    output = []
    with open(file, 'r') as handle:
        assembly = SeqIO.read(handle, 'fasta')
        for key in block_dict.keys():
            if key.find(prefix) > -1:
                for block in block_dict[key].keys():
                    start, stop, strand = block_dict[key][block]
                    seq = assembly.seq if strand == '+' else assembly.seq.reverse_complement()
                    block_seq = str(seq)[int(start) - 1: int(stop)]
                    record = SeqRecord(Seq(block_seq), id=prefix, name=block,  description=f'syntenic block {block}')
                    output.append(record)
    return output



@click.command()
@click.option('--input', '-i', help='Specify an input file', type=click.Path(exists=True))
@click.option('--assemblies', '-a', help='Specify a directory containing target assemblies', type=click.Path(exists=True))
@click.option('--output', '-o', help='Specify an output directory', type=click.Path(exists=False))

def execute(input, assemblies, output):
    with open(input, 'r') as handle:
            click.echo('The provided file is a PaReBrick sorted block list')
            click.echo('Parsing PaReBrick output...')
            block_dict = parseParebrickBlocks(handle)
            click.echo(len(block_dict.keys()))
    if not os.path.exists(output):
        os.makedirs(output)
    for file in os.listdir(assemblies):
        prefix = file.split('.')[0]
        file = os.path.join(assemblies, file)
        block_seqs = extractBlocks(file, block_dict)
        output_dir = os.path.join(output, prefix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            for block in block_seqs:
                output_file = os.path.join(output_dir, f'{block.name}.fna')
                with open(output_file, 'w') as handle:
                    SeqIO.write(block, handle, 'fasta')


if __name__ == '__main__':
    execute()
