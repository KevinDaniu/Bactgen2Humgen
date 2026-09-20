#!/bin/bash

for strain in $(ls ../result/bactgene2humgene/); do
    echo "========== 开始 $strain =========="
    bash 3.2-run_strain.sh "$strain"
    echo "========== 结束 $strain =========="
done
