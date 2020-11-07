#!/usr/bin/env python3

from Bio import Align
from Bio import AlignIO
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from distMat import openAlignment

import click
import os

"""
Concatenates alignments containing the same number of entries and shared sequence identifiers. Returns the concatenated alignment and a list of coordinates.
"""


##### consider changing it
#def dirwise(dir):
#    """Applies function to all the files in the directory"""
#    def apply_dirwise(func):
#        def wrapper(*args, **kwargs):
#            output_list = []
#            files = os.files(dir)
#            for file in alignments:
#                output_list.append(func(file, *args, **kwargs))
#            return output_list
#        return wrapper

##### consider changing it


def openAlignmentWithName(file:str, format:str)->tuple:
    aln = openAlignment(file, format)
    prefix='.'.join(file.split('.')[:-1])
    return (prefix, aln)


def openAlignmentDirwise(dir:str, format:str)->list:
    """
    Applies 'openAlignmentWithNames' function to all the files in the directory.
    Consider rewriting it wth decorators the other day
    """
    output_list = []
    files = os.listdir(dir)
    for file in files:
        output_list.append(openAlignmentWithName(os.path.join(dir, file), format))
    return output_list

def concatenateAlignments(alignments:list)->tuple:
    concatenate_dict = defaultdict(list)
    metadata = []
    start = 1
    for name, aln in alignments:
        aln.sort()
        length = aln.get_alignment_length()
        name = name.split('/')[-1].split('.')[0]
        for record in aln:
            concatenate_dict[record.id.split('_')[0]].append(str(record.seq))
        metadata.append(f'DNA, {name} = {start}-{start+length-1}')
        start = start + length
    return (Align.MultipleSeqAlignment(SeqRecord(Seq(''.join(seq), alphabet=SingleLetterAlphabet()), id=id) for (id,seq) in concatenate_dict.items()), metadata)


@click.command()
@click.option('--input', '-i', help='Specify the path to the input directory', type=click.Path(exists=True))
@click.option('--format', '-f', help='Specify an alignment format (default is "fasta")', default='fasta')
@click.option('--output_alignment', '-oa', help='Specify the path to write the concatenated alignment to', type=click.Path(exists=False))
@click.option('--output_metadata', '-om', help='Specify the path to the alignment layout', type=click.Path(exists=False), default=None)

def execute(input, format, output_alignment, output_metadata):
    print('Initializing...')
    alignments = openAlignmentDirwise(input, format)
    concatenated = concatenateAlignments(alignments)
    with open(output_alignment, 'w') as handle:
        AlignIO.write(concatenated[0], handle, 'fasta')
    if output_metadata:
        with open(output_metadata, 'w') as handle:
            handle.write('\n'.join(concatenated[1]))

if __name__ == '__main__':
    execute()
