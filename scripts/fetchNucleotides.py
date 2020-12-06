#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict

import click
import gzip
import os
import pandas as pd

"""
Fetches the nucleotide sequences from the assembly annotations by the list of accessions from an another file (e.g., protein FASTA file)
For the sake of compatibility with different annotation formats, inclusion check is performed instead of name identity check
"""

def parseOrthogroups(input:str)->dict:
    """
    Parses OrthoFinder orthogroup table and returns a dictionary containing accessions for each
    single-copy ortholog in each genome
    """

    out = defaultdict(dict)
    with open(input, 'r') as handle:
        og = pd.read_csv(handle, sep='\t', header=0)
        genomes = list(filter(lambda x: x.find('Orthogroup') == -1, og.columns.values))
        for i in range(0, og.shape[0]):
            orths = pd.Series(og.iloc[i][1:])
            orths = orths[-orths.isnull()]
            orths = list(filter(lambda x: x.find('Orthogroup') == -1, orths))
            if len(orths) == len(genomes):
                for genome in genomes:
                    out[genome][og[genome][i]] = og.iloc[i, 0]

    return out


def findFile(directory:str, *args):
    output_list = []
    for path, dirs, files in os.walk(directory):
        for file in files:
            if all(x in file for x in args):
                output_list.append(os.path.join(path, file))
        for dir in dirs:
            findFile(dir, *args)
    return output_list

def extractOrthologousGenes(indict:defaultdict, input:str,  output:str, mode:str, misc:list=[]):
    """
    Extract single-copy orthologs according to the parseOrthogroups output and saves them
    either filewise or genomewise. NOTE: the implementation of findFile here assumes that 
    a provided set of identifiers will unequivocally identify s CDS-containing file
    """

    if not os.path.exists(output):
        os.makedirs(output)

    for gen in indict.keys():
        pattern = '_'.join(gen.split('_')[0:2])
        record_list = []
        file = findFile(input, *[pattern, *misc])
        with gzip.open(*file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                for key in indict[gen].keys():
                    if key in record.id:
                        record.id = key
                        if mode == 'orth':
                            if not os.path.exists(os.path.join(output, pattern)):
                                os.makedirs(os.path.join(output, pattern))
                            record.description = f'genome={pattern} | orthogroup={indict[gen][key]}'
                            with open(os.path.join(output, pattern, f'{indict[gen][key]}_nuc.fna'), 'w') as handle:
                                SeqIO.write([record], handle, 'fasta')
                        elif mode == 'genome':
                            record.description = f'orthogroup={indict[gen][key]}'
                            record_list.append(record)
        if mode == 'genome':
            with open(os.path.join(output, f'{pattern}_nuc.fna'), 'w') as handle:
                SeqIO.write(record_list, handle, 'fasta')


@click.command()
@click.option('--input', '-i', help='A path to the directory containing nucleic acid sequences', type=click.Path(exists=True))
@click.option('--table', '-t', help='A path to an OrthoFinder ("Orthogroups.tsv") orthogroup table', type=click.Path(exists=True))
@click.option('--options', '-p', help='Search terms used to identify the corresponding NA file', multiple=True, type=str, default='')
@click.option('--mode', '-m', help='Define whether the sequences should be extracted in a genome-wise (mode="genome") or ortholog-wise (mode="orth") manner. In the latter case, a directory named after the genome accession is created containing a single file for each orthogroup', type=click.Choice(['genome', 'orth'], case_sensitive=False), default='genome')
@click.option('--output', '-o', help='A path to the output directory', type=click.Path(exists=False))

def execute(input, output, table, mode,  options):
    og = parseOrthogroups(table)
    extractOrthologousGenes(og, input, output, mode, options)


if __name__ == '__main__':
    execute()
