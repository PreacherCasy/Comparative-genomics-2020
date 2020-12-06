#!/usr/bin/env python3

import click
import os
import pandas as pd

"""
"""

def replaceLabels(input, mapper):
    newlines = []
    for line in input.readlines():
        for i in range(mapper.shape[0]):
            line = line.replace(mapper.iloc[i,1], mapper.iloc[i, 0])
        newlines.append(line)
    return newlines

@click.command()
@click.option('--input', '-i', help='', type=click.Path(exists=True))
@click.option('--mapper', '-m', help='', type=click.Path(exists=True))
@click.option('--output', '-o', help='', type=click.Path(exists=False))

def execute(input, mapper, output):
    with open(input, 'r') as filehandle:
        with open(mapper, 'r') as maphandle:
            map_df = pd.read_csv(maphandle, sep='\t', header=None)
            out_string = replaceLabels(filehandle, map_df)
    with open(output, 'w') as handle:
        for line in out_string:
            handle.write(line)


if __name__ == '__main__':
    execute()

