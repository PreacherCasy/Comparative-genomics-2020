#!/usr/bin/env python3

from extractModel import extractModel

import click
import _io
import os

"""
Prepares a partition file by extracting model from the modeltest-ng and adding them to the partioning scheme
"""

def extractAll(file:_io.TextIOWrapper, pred:str)->list:
    output_list = []
    while file:
        try:
            record = extractModel(file, pred)
            output_list.append(record)
        except:
            break
    return output_list


@click.command()
@click.option('--input', '-i',
 help='Specify an input file ("*.out" from the modeltest-ng output)', type=click.Path(exists=True))
@click.option('--criterion', '-c', help='Select a criterion to choose models by. Available options are AIC, AICc and BIC (default)',
 type=click.Choice(['AIC', 'AICc', 'BIC']), default='BIC')
@click.option('--output', '-o', help='Specify the path to the input file',
 type=click.Path(exists=False), default=None)

def execute(input, criterion, output):
    with open(input, 'r') as handle:
        output_list = extractAll(handle, criterion)
    if output:
        with open(output, 'w') as handle:
            handle.write('\n'.join(output_list))
    else:
        click.echo('\n'.join(list(map(lambda x: 
                  f'The best model according to {criterion} is {x}', output_list))), err=False)

if __name__ == '__main__':
    execute()
