#! /usr/bin/env bash
# Auth: Jennifer Chang
# Date: 2020/02/03

set -e
set -u

CODEDIR=~/Desktop/focus/Swine_Survey/code

ARR=(HH NN PB2 PB1 PA NP M NS)

for SEG in ${ARR[@]}
do
    echo $SEG
    ${CODEDIR}/batchFetchGB.sh ${SEG}.ids > ${SEG}.gb
    ./gb2seqNT.pl ${SEG}.gb > ${SEG}.fna
done