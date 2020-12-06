#!/usr/bin/env python3

from Bio import SeqIO
from collections import ChainMap
from collections import defaultdict
from typing import Generator
from fetchNucleotides import findFile

import click
import gzip
import os

"""
Sorts prokaryotic genomes into chromosomes, plasmids and unattributed contigs/fragments
and writes them respective to their type
"""

def fetchContigTypes(file:Generator, alias:str, allowed_types:list, keep:bool)->defaultdict:
    """
    Sorts all contigs found according to the listed types.
    If keep is is set to 'True', assigns all the contigs failing to match the specified types to the
    'misc' category
    """

    output = defaultdict(lambda: defaultdict(list))
    for record in file:
        descr = record.description
        if sum(list(map(lambda x: max(descr.find(x),
         descr.find(x.capitalize())), allowed_types))) > -1*len(allowed_types):
            for allowed_type in allowed_types:
                if descr.find(allowed_type) > -1:
                    output[allowed_type][alias].append(record)
        elif keep:
            output['misc'][alias].append(record)
    return output


@click.command()
@click.option('--input' , '-i', help='Specify the path to the input file',
 type=click.Path(exists=True))
@click.option('--allowed-types', '-t', help='Specify the contig types to be sought for',
type=list, default=['chromosome', 'contig', 'plasmid', 'scaffold'])
@click.option('--keep', '-k', help='''
Specify whether the contigs failing to meet any name should be considered (default is False).
If enabled, the unnamed contigs are stored under the "misc" name.
''', is_flag=True)
@click.option('--output', '-o', help='Specify the path to the output directory', type=click.Path(exists=False))

def execute(input, allowed_types, keep, output):
    files = [x for x in findFile(input, 'GCF', 'genomic', 'fna') if x.find('cds') == -1]
    for file in files:
        alias = file.split('/')[-1].split('.')[0]
        with gzip.open(file, 'rt') as handle:
            assembly = SeqIO.parse(handle, 'fasta')
            assembly_contigs = fetchContigTypes(assembly, alias, allowed_types, keep)
        if not os.path.exists(output):
            os.makedirs(output)
        for key in assembly_contigs.keys():
            contigwise = os.path.join(output, key)
            if not os.path.exists(contigwise):
                os.makedirs(contigwise)
            for assembly in assembly_contigs[key].keys():
                file = os.path.join(contigwise, assembly)
                with open(f'{file}.fna', 'w') as handle:
                    SeqIO.write(assembly_contigs[key][assembly], handle, 'fasta')



if __name__ == '__main__':
    execute()
