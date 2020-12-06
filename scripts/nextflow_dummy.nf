#!/usr/bin/env nextflow

params.indir = '/home/yura/compgenomics_assignment'


process listDir {

output:
stdout result

"""
ls ${params.indir}
"""

}

result.view { it }
