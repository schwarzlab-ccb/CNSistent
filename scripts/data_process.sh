#!/bin/bash
threads=30

data="../data"
out="../out"
assembly="hg19"

set -x

cd "$(dirname "$0")" # Set path to the script's path

mkdir -p $out

segment whole --remove gaps --out ${out}/gaps_hg19_segs.bed --verbose --assembly hg19
segment whole --remove gaps --out ${out}/gaps_hg38_segs.bed --verbose --assembly hg38

for dataset in TRACERx PCAWG TCGA_hg19;
do
    echo "Processing $dataset with assembly $assembly"      
    common_args="--threads $threads --verbose --assembly $assembly"
    cns align "${data}/${dataset}_cns_raw.tsv" --samples "${data}/${dataset}_samples_raw.tsv" --out "${out}/${dataset}_cns_align.tsv" $common_args
    cns infer "${out}/${dataset}_cns_align.tsv" --samples "${data}/${dataset}_samples_raw.tsv" --out "${out}/${dataset}_cns_imp.tsv" $common_args
    cns coverage "${out}/${dataset}_cns_align.tsv" --samples "${data}/${dataset}_samples_raw.tsv" --out "${out}/${dataset}_samples_align.tsv" $common_args        
    cns coverage "${out}/${dataset}_cns_align.tsv" --samples "${data}/${dataset}_samples_raw.tsv" --out "${out}/${dataset}_samples.tsv" $common_args --segments "${out}/gaps_${assembly}_segs.bed"
    cns ploidy "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${out}/${dataset}_samples.tsv" $common_args --segments "${out}/gaps_${assembly}_segs.bed"
    cns breakage "${out}/${dataset}_cns_imp.tsv" --samples "${out}/${dataset}_samples.tsv" --out "${out}/${dataset}_samples.tsv" $common_args --segments "${out}/gaps_${assembly}_segs.bed"
done

