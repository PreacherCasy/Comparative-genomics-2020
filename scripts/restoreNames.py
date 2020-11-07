#!/usr/bin/env python3

import click
import os
import pandas as pd

"""
Restores genome-to-sequence accordance based on the Roary presence/absence matrix
"""


def openDF(input:str, sep:str='\t', header:int=0):
    with open(input, 'r') as handle:
        og = pd.read_csv(handle, sep=sep, header=header)
    return og


def namesFromColumns(df:pd.core.frame.DataFrame, ncol:tuple=(0,0), pattern:str=None)->dict:
    output_dict = dict()
    if pattern:
        genomes = list(filter(lambda x: x.find(pattern) > -1, df.columns.values))
    else:
        genomes = df.columns.values[slice(*ncol)]
#    gen_ind = [df.columns.values.index(x) for x in genomes]
    for genome in genomes:
        output_dict[genome] = df[genome][0].split('.')[0].split('_')[0]
    return output_dict

@click.command()
@click.option('--input', '-i', help='Specify an input orthogroup table', type=click.Path(exists=True))
@click.option('--sep', '-sp', help='Table field separator', type=str, default='\t')
@click.option('--header', '-h', help='Table header specification', type=int, default=None)
@click.option('--ncol', '-nc', help='Columns to be used as genome specifications', type=tuple, default=(0,0))
@click.option('--pattern', '-p', help='Pattern to find genome-identifying columns by. Overrides "ncol" argument', type=str, default=None)

def execute(input, sep, header, ncol, pattern):
    df = openDF(input, sep, header)
    name_dict = namesFromColumns(df, ncol, pattern)
    print(name_dict)

if __name__ == '__main__':
    execute()
