#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.segseqs="../data/*.fna"
params.outdir="./Results"

process mafft {
    publishDir "${params.outdir}/01_Align"

    input: path(segment_fasta)
    output: path("*_aln.fna")

    script:
    """
    #! /usr/bin/env bash
    mafft --auto ${segment_fasta} > ${segment_fasta.simpleName}_aln.fna
    """
}

process fasttree {
    publishDir "${params.outdir}/02_Trees"
    input: path(aln_fna)
    output: path("*.tre")
    script:
    """
    #! /usr/bin/env bash
    FastTreeMP -nt $aln_fna > ${aln_fna.simpleName}.tre
    """
}

workflow {
    seg_ch = channel.fromPath(params.segseqs)
    seg_ch | mafft | fasttree
}