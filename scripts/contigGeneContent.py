#!/usr/bin/env python3

from BCBio import GFF
from collections import defaultdict

import click
import gzip
import os
import pandas as pd

"""
"""

def contigGeneContent(input:click.Path)->pd.core.frame.DataFrame:
    assembly = input.split('/')[-1].split('.')[0]
    output_dict = defaultdict(list)
    with gzip.open(input, 'rt') as handle:
        gff = GFF.parse(handle, limit_info={'gff_type': ['region', 'gene']})
        for rec in gff:
            try:
                output_dict[rec.id].append(rec.features[0].qualifiers['genome'][0])
            except KeyError:
                output_dict[rec.id].append('NA')
            output_dict[rec.id].append(0)
            for feature in rec.features:
                if feature.type == 'gene' and 'protein_coding' in feature.qualifiers['gene_biotype']:
                    output_dict[rec.id][1] += 1
    output = pd.DataFrame({'Assembly': [assembly]*len(output_dict.keys()),
                           'Contig': list(output_dict.keys()),
                           'ContigType': list(map(lambda x: x[0], output_dict.values())),
                           'NoOfGenes': list(map(lambda x: x[1], output_dict.values()))})
    return output


@click.command()
@click.option('--input', '-i', help='Specify a path to an input GFF file', 
type=click.Path(exists=True))
@click.option('--output', '-o', help='Specify a path to write the output table to',
type=click.Path(exists=False))

def execute(input, output):
    df = contigGeneContent(input)
    with open(output, 'w') as handle:
        df.to_csv(handle, index=False, header=True)


if __name__ == '__main__':
    execute()
