#!/usr/bin/env python3

from alignClusters import alignFilewise

import click
import os
import subprocess

"""
Aligns all sequence files given in a directory and writes output to the specified directory
"""

@click.command()
@click.option('--input', '-i', help='An input directory containing alignment files', type=click.Path(exists=True))
@click.option('--output', '-o', help='An output directory', type=click.Path(exists=False))
@click.option('--aligner', '-a', help='Choose between PRANK or MAFFT to perform MSA. Default is MAFFT. Note that the aligner should be reachable from $PATH', type=click.Choice(['mafft', 'prank'], case_sensitive=True), default='mafft')
@click.option('--options', '-m', help='Options for the aligner chosen', type=str, multiple=True)

def execute(input:str, output:str, aligner, options):
    if not os.path.exists(output):
        os.makedirs(output)
    print(f'Options are: {options}')
    files = os.listdir(input)
    for file in files:
        try:
            alignFilewise(os.path.join(input, file), output, aligner, *options)
        except subprocess.CalledProcessError as e:
            click.echo(f'{e}: {os.path.join(input, file)} is not a multifasta file!')
            continue

if __name__ == '__main__':
    execute()
