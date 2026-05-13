#!/bin/bash

set -x

# set current directory to the location of the script
cd "$(dirname "$0")"

threads=10
subsplit=5
data="../data"
out="../out"

segment whole --out "${out}/segs_20MB_full.bed" --split 20000000
segment whole --out "${out}/segs_20MB.bed" --split 20000000 --remove gaps --filter 2000000 
segment whole --out "${out}/segs_10MB.bed" --split 10000000 --remove gaps --filter 1000000 
segment whole --out "${out}/segs_5MB.bed" --split 5000000 --remove gaps --filter 500000
segment whole --out "${out}/segs_3MB.bed" --split 3000000 --remove gaps --filter 300000  
segment whole --out "${out}/segs_2MB.bed" --split 2000000 --remove gaps --filter 200000
segment whole --out "${out}/segs_1MB.bed" --split 1000000 --remove gaps --filter 100000
segment whole --out "${out}/segs_500KB.bed" --split 500000 --remove gaps --filter 50000
segment whole --out "${out}/segs_250KB.bed" --split 250000 --remove gaps --filter 25000
segment whole --out "${out}/segs_100KB.bed" --split 100000 --remove gaps --filter 10000
segment "${data}/COSMIC_consensus_genes.bed" --out "${out}/segs_COSMIC.bed" 
segment "${data}/ENSEMBL_coding_genes.bed" --out "${out}/segs_ENSEMBL.bed"
segment whole --out "${out}/segs_whole.bed"
segment "arms" --out "${out}/segs_arms_full.bed"
segment "arms" --out "${out}/segs_arms.bed" --remove gaps --filter 1000000
segment "bands" --out "${out}/segs_bands.bed" --remove gaps --filter 100000
for dist in 1MB 500KB 250KB; do
    python ./data_cluster.py $dist
done


# TRACERx PCAWG TCGA_hg19
for dataset in TRACERx PCAWG TCGA_hg19; 
do    
    echo "Processing $dataset"
    shared_args="${out}/${dataset}_cns_imp.tsv --samples ${out}/${dataset}_samples.tsv  --verbose --threads $threads --subsplit $subsplit"
    cns aggregate --segments "${out}/segs_20MB.bed" --out "${out}/${dataset}_bin_20MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_20MB_full.bed" --out "${out}/${dataset}_bin_20MB_full.tsv" $shared_args    
    cns aggregate --segments "${out}/segs_whole.bed" --out "${out}/${dataset}_bin_whole.tsv" $shared_args
    cns aggregate --segments "${out}/segs_arms.bed" --out "${out}/${dataset}_bin_arms.tsv" $shared_args
    cns aggregate --segments "${out}/segs_arms_full.bed" --out "${out}/${dataset}_bin_arms_full.tsv" $shared_args
    cns aggregate --segments "${out}/segs_bands.bed" --out "${out}/${dataset}_bin_bands.tsv" $shared_args
    cns aggregate --segments "${out}/segs_20MB.bed" --out "${out}/${dataset}_bin_20MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_10MB.bed" --out "${out}/${dataset}_bin_10MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_5MB.bed" --out "${out}/${dataset}_bin_5MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_3MB.bed" --out "${out}/${dataset}_bin_3MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_2MB.bed" --out "${out}/${dataset}_bin_2MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_1MB.bed" --out "${out}/${dataset}_bin_1MB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_500KB.bed" --out "${out}/${dataset}_bin_500KB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_250KB.bed" --out "${out}/${dataset}_bin_250KB.tsv" $shared_args
    cns aggregate --segments "${out}/segs_COSMIC.bed" --out "${out}/${dataset}_bin_COSMIC.tsv" $shared_args
    # Uncomment for finer segmentations, but beware of long runtimes and large output files
    # cns aggregate --segments "${out}/segs_100KB.bed" --out "${out}/${dataset}_bin_100KB.tsv" $shared_args
    # cns aggregate --segments "${out}/segs_ENSEMBL.bed" --out "${out}/${dataset}_bin_ENSEMBL.tsv" $shared_args
    for dist in 1MB 500KB 250KB; do
        cns aggregate --segments "${out}/segs_merge_${dist}.bed" --out "${out}/${dataset}_bin_merge_${dist}.tsv" $shared_args 
    done
    # Uncomment for comparison with other aggregation methods
    # cns aggregate --segments "${out}/segs_COSMIC.bed" --out "${out}/${dataset}_bin_COSMIC_min.tsv" --how min $shared_args
    # cns aggregate --segments "${out}/segs_COSMIC.bed" --out "${out}/${dataset}_bin_COSMIC_max.tsv" --how max $shared_args
done
