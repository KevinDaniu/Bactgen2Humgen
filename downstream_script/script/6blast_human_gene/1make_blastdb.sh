#!/bin/bash

mkdir -p result/blastdb_human
makeblastdb -in rawdata/gencode.v49.pc_transcripts.fa -dbtype nucl -out result/blastdb_human/transcript_db
