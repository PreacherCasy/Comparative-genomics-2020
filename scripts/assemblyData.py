#!/usr/bin/env python3

from Bio import Entrez
from Bio import SeqIO
from collections import defaultdict
from download_assemblies import get_summary

import click
import os
import pandas as pd

@click.command()
@click.option('--email', '-e', help='Your email to introduce to Entrez', type=str)
@click.option('--name', '-n', help='Query organism name', type=str)
@click.option('--output', '-o', help='Output directory name', type=click.Path(exists=False), default=None)


def writeStatistics(email:str, name:str, output:str):

    """
    Writes assembly status for each entry found by the query
    """

    Entrez.email, Entrez.tool = email, 'MyCustomScript'

    search_output = Entrez.read(Entrez.esearch(db='assembly', term=name, retmax=100000))['IdList']
    prefix = f'{"_".join(name.split(" "))}'
    full_dict = defaultdict(list)
    summary_dict = dict()


    id_list = ','.join(search_output)
    summaries = get_summary(id_list)
    for summary in summaries['DocumentSummarySet']['DocumentSummary']:
        url = summary['FtpPath_RefSeq']
        if url != '':
           full_dict['Accession'].append(summary['AssemblyAccession'])
           full_dict['Name'].append(summary['AssemblyName'])
           full_dict['Level'].append(summary['AssemblyStatus'])
           full_dict['Coverage,x'].append(summary['Coverage'])

    out_full = pd.DataFrame(full_dict)
    out_full.to_csv(os.path.join(output, f'{prefix}_full_report.tsv'), sep='\t', index=False, header=True)

if __name__ == '__main__':
    writeStatistics()
