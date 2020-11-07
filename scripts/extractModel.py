#!/usr/bin/env python3

import click
import os

"""
Extracts best model according to modeltest-ng
"""

def extractModel(file:str, pred:str)->str:
#    model = ''
    for line in file:
        if line.find(f'Best model according to {pred}') > -1:
            next(file)
            model = next(file).rstrip('\n')
            model = model.replace(' ', '').split(':')[1]
            break
    return model


@click.command()
@click.option('--input', '-i', help='Specify a modeltest-ng "out" file to process', type=click.Path(exists=True))
@click.option('--criterion', '-c', help='Select a criterion to choose a model by. Available options are AIC, AICc and BIC (default)', type=click.Choice(['AIC', 'AICc', 'BIC']), default='BIC')

def execute(input, criterion):
    with open(input, 'r') as handle:
        print(f'The best model according to {criterion} is {extractModel(handle, criterion)}')

if __name__ == '__main__':
    execute()
