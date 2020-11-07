#!/usr/bin/env python3

import click
import os

"""
Parse ModelTest-NG results and retrieves the model specification suitable for further upload to RAxML-NG 
"""

def parseModelTestNg(input:str)->str:
    with open(input, 'r') as handle:
        all_lines = handle.readlines()
        print(len(all_lines))
        print(all_lines.index('Selection options:\n'))


@click.command()
@click.option('--input', '-i', help='A path to the input file', type=click.Path(exists=True))

def execute(input):
    parseModelTestNg(input)

if __name__ == '__main__':
    execute()
