#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
from fetchNucleotides import findFile, parseOrthogroups

import click
import gzip
import math
import os
import pandas as pd
import zipfile

"""
Fetches nucleotide sequences for a specified OrthoFinder cluster
"""

def dictByQuery(og:pd.core.frame.DataFrame, ncols:tuple=(1,1), *args)->defaultdict:
    output_dict = defaultdict(dict)
    for query in args:
        for cluster in query.split(','):
#            print(cluster)
            idx = pd.Index(og.iloc[:,0]).get_loc(cluster)
            ncols = tuple(map(lambda x: int(x), ncols))
            for i in range(*ncols):
                if isinstance(og.iloc[idx,i], str):
                    output_dict[cluster][og.columns.values[i]] = og.iloc[idx,i]
#    print(output_dict.keys())
    return output_dict


def fetchByQuery(input:defaultdict, dir:str, misc:tuple={}):
    record_dict = defaultdict(list)
    for key1 in input.keys():
#        print(key1)
        for key2 in input[key1]:
#            print(key2)
            file = findFile(dir, *(key2, *misc))
            with gzip.open(file, 'rt') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    if record.id.find(input[key1][key2]) > -1:
                        record.id = key2
                        record_dict[key1].append(record)
#        print(record_dict.keys())
    return record_dict


@click.command()
@click.option('--table', '-t', help='', type=click.Path(exists=True))
@click.option('--input', '-i', help='', type=click.Path(exists=True))
@click.option('--query', '-q', help='', type=str, multiple=True)
@click.option('--output', '-o', help='', type=click.Path(exists=False), default='.')
@click.option('--ncols', '-nc', help='', type=(int, int), default=(1,1))
@click.option('--pattern', '-p', help='', type=str, multiple=True, default='')

def execute(input, table, output, ncols, query, pattern):
    if not os.path.exists(output):
        os.makedirs(output)
    with open(table, 'r') as handle:
       df = pd.read_csv(table, header=0, sep='\t')
    ogs = dictByQuery(df, ncols, *query)
    output_dict = fetchByQuery(ogs, input, pattern)
    print(output_dict.keys())
    for key in output_dict.keys():
        SeqIO.write(output_dict[key], os.path.join(output, f'{key}.fna'), 'fasta')

if __name__ == '__main__':
    execute()
