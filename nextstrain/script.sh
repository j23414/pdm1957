#! /usr/bin/env bash

set -eu
#>CY021709|A/AA/Huston/1945|H1N1|Human|USA|1945

ARR=(HH NN PB1 PB2 PA NP M NS)
DATA_DIR=../data
REFERENCE_NAME="A/Brevig_Mission/1/1918"

mkdir -p results
for seg in "${ARR[@]}"; do
    echo "Processing segment $seg"
    mkdir -p results/$seg
    echo "accession|strain|subtype|host|country|date" \
    | tr '|' '\t' \
    > results/$seg/metadata.tsv
    cat $DATA_DIR/${seg}_annot.fna \
    | smof grep "|" \
    | grep ">" \
    | sed 's/>//g' \
    | tr '|' '\t' \
    >> results/$seg/metadata.tsv

    tsv-uniq -H -f strain results/$seg/metadata.tsv > temp.txt
    mv temp.txt results/$seg/metadata.tsv

    cat $DATA_DIR/${seg}_annot.fna \
    | smof grep "|" \
    | sed 's/>[^|]*|/>/g' \
    | sed 's/|.*//g' \
    | smof uniq --first-header \
    > results/$seg/sequences.fasta

    if [[ ! -f "results/$seg/aligned.fasta" ]]; then
        echo "Aligning segment $seg"
        augur align \
        --sequences results/$seg/sequences.fasta \
        --reference-name ${REFERENCE_NAME} \
        --output results/$seg/aligned.fasta
    fi

    if [[ ! -f "results/$seg/tree-raw.nwk" ]]; then
        echo "Building tree for segment $seg"
        augur tree \
        --alignment results/$seg/aligned.fasta \
        --output results/$seg/tree-raw.nwk
    fi

    if [[ ! -f "results/$seg/tree.nwk" ]]; then
       echo "Refine tree for segment $seg"
       augur refine \
        --tree results/$seg/tree-raw.nwk \
        --alignment results/$seg/aligned.fasta \
        --metadata results/$seg/metadata.tsv \
        --output-tree results/$seg/tree.nwk \
        --output-node-data results/$seg/branch_lengths.json \
        --timetree \
        --coalescent opt \
        --date-confidence \
        --date-inference marginal \
        --root ${REFERENCE_NAME}
    fi

    if [[ ! -f "results/$seg/traits.json" ]]; then
      augur traits \
        --tree results/$seg/tree.nwk \
        --metadata results/$seg/metadata.tsv \
        --output results/$seg/traits.json \
        --columns accession subtype host country date \
        --confidence
    fi

    if [[ ! -f "results/$seg/nt_muts.json" ]]; then
      augur ancestral \
        --tree results/$seg/tree.nwk \
        --alignment results/$seg/aligned.fasta \
        --output-node-data results/$seg/nt_muts.json \
        --inference joint
    fi

    if [[ ! -f "results/$seg/aa_muts.json" ]]; then
      augur translate \
        --tree results/$seg/tree.nwk \
        --ancestral-sequences results/$seg/nt_muts.json \
        --reference-sequence $DATA_DIR/${seg}_annot.fna \
        --output-node-data results/$seg/aa_muts.json
    fi

    augur export v2 \
        --tree results/$seg/tree.nwk \
        --metadata results/$seg/metadata.tsv \
        --node-data \
          results/$seg/branch_lengths.json \
          results/$seg/traits.json \
          results/$seg/nt_muts.json \
          results/$seg/aa_muts.json \
        --color-by-metadata subtype host country \
        --output auspice/$seg.json

done