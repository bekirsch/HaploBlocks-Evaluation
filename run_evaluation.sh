#!/bin/bash

## threads=$1 is an inputparameter specifying the number of threads to run

# Check how many cores are physically available
avail=$(nproc --all)
# Prevent using more threads than available cores
if (( $1 > $avail )); then
  threads=$(($avail-1))
else
  threads=$1
fi

if [[ -z "$*" ]]; then
 threads=1
fi

echo 'Running script with' $threads 'threads'

# Make directory for simulations
mkdir results/evaluation/Additive_10Mb_10kNe

# Define function for simulating
simulating() {
# Setting random seed
seed=$(openssl rand 4 | od -DAn);

# Running SLiM
slim -s $seed scripts/Additive.slim &>/dev/null;

trees=$(echo "results/evaluation/Additive_10Mb_10kNe/simulation${seed}/${seed}_sC0.02_mF*.trees" | tr -d ' ');

for file in $trees; do
  python3 scripts/recapitation.py -i $file &>/dev/null;
  line=$(cat ${file/.trees/.trees.vcf} | grep -n '4999999' | cut -f1 | cut -d":" -f1)
  cat ${file/.trees/.trees.vcf} | awk -F '\t' -v OFS='\t' -v m=$line -v n=4 -v el='0' 'NR == m { $n = el } 1' | awk -F '\t' -v OFS='\t' -v m=$line -v n=5 -v el='1' 'NR == m { $n = el } 1' | gzip > ${file/.trees/.trees.uniform.vcf.gz}
  rm ${file/.trees/.trees.vcf}
done
}
export -f simulating

# Run simulations in multiple threads
parallel -j $threads --no-notice simulating ::: {1..50}

# creating lookup-table
tools/haploblocks/filter_lookup -max_k 2000 > results/evaluation/Additive_10Mb_10kNe/ancestry.lookup

# Make directory for output
mkdir results/evaluation/Additive_10Mb_10kNe/output

# Define function for running haploblocks
haploblocks() {
vcf_gz=$1;
vcf=${vcf_gz/.vcf.gz/.vcf}

zcat $vcf_gz > $vcf

cmap=${vcf/.vcf/.vcf.positions};
rmap=${cmap/.positions/.positions.map};

### Extract positions and generate genetic map ###
tools/haploblocks/extract_positions -i $vcf -o $cmap &>/dev/null;

awk -v OFS='\t' '{print "chr1", "snp"NR, (50*log(1/(1-(2*1e-8*$0)))), $0}' $cmap | tr ',' '.' > $rmap;

tools/haploblocks/full --out_folder results/evaluation/Additive_10Mb_10kNe/output --vcf_path $vcf --genetic_map_path $rmap --lookup_path results/evaluation/Additive_10Mb_10kNe/ancestry.lookup --remove &>/dev/null;

### Clean Up ###
rm $vcf;
rm $cmap;
rm $rmap;
}
export -f haploblocks

# Run haploblocks in multiple threads
find results/evaluation/Additive_10Mb_10kNe/simulation*/*.uniform.vcf.gz | parallel -j $threads --no-notice haploblocks {}

# Count simulations
all=$(ls results/evaluation/Additive_10Mb_10kNe/output/*filtered.sHat.csv | wc -l)
files=$((all / 13))

# Plot results
Rscript scripts/Plot_evaluation.R results/evaluation/Additive_10Mb_10kNe/output $files
