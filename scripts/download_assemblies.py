#!/usr/bin/env python3

import click
import os
import urllib
from Bio import Entrez
from Bio import SeqIO

def get_summary(name):

    esummary_handle = Entrez.esummary(db="assembly", id=name, report="full", retmode='xml')
    esummary_record = Entrez.read(esummary_handle, validate=False)

    return esummary_record


@click.command()
@click.option('--email', '-e', help='Your email to introduce to Entrez', type=str)
@click.option('--name', '-n', help='Query organism name', type=str)
@click.option('--output', '-o', help='Output directory name', type=click.Path(exists=False))


def downloadAssemblies(email:str, name:str, output:str):

    Entrez.email, Entrez.tool = email, 'MyCustomScript'

    if not os.path.exists(output):
        os.makedirs(output)

    search_output = Entrez.read(Entrez.esearch(db='assembly', term=name, retmax=100000))['IdList']

    id_string = ','.join(search_output)
    dir_names = {'Complete Genome':  'complete_assemblies', 'Chromosome': 'chromosome_level'}
    summaries = get_summary(id_string)
    for summary in summaries['DocumentSummarySet']['DocumentSummary']:
        url = summary['FtpPath_RefSeq']
        if url != '':
            label=os.path.basename(url)
            files = ['_genomic.fna.gz', '_cds_from_genomic.fna.gz', '_protein.faa.gz', '_genomic.gff.gz',  '_genomic.gbff.gz']
            level = summary['AssemblyStatus']
            if level in ['Complete Genome', 'Chromosome']:
                new_output = os.path.join(output, dir_names[level])
                if not os.path.exists(new_output):
                    os.makedirs(new_output)
                os.makedirs(os.path.join(new_output, label))
                for file in files:
                    link = os.path.join(url,f'{label}{file}')
                    filepath = os.path.join(new_output, label)
                    urllib.request.urlretrieve(link, os.path.join(filepath,  f'{label}{file}'))


if __name__ == '__main__':
    downloadAssemblies()
