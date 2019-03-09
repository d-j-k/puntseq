#!/usr/bin/env sh

# RUN FROM ANALYSIS DIRECTORY

mkdir -p data/loman/basecalled

URL=https://nanopore.s3.climb.ac.uk/Zymo-GridION-EVEN-BB-SN.fq.gz

wget "$URL" -O data/loman/basecalled/loman_all_passed.fastq.gz
