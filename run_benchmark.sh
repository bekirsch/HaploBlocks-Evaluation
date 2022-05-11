#!/bin/bash

cd chr/

# If you change this filename, make shure to provide the vcf.gz OR skip the next command (zcat ${file/.vcf/.vcf.gz} > $file)!
# Also genetic map must be provided in the same folder and needs to be changed accordingly (line 72)
file1=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
file=chr22.1kgpp3.vcf

# Unzip the VCF.gz
zcat $file1 > $file

# Remove indels and filter for MAF = 0.01
vcftools --vcf $file --remove-indels --maf 0.01 --recode --out ${file/.vcf/.vcf.MAF0.01} &>/dev/null

# Remove the VCF again
rm $file

# Extract a list of SNPs
cat ${file/.vcf/.vcf.MAF0.01.recode.vcf} | grep -v '#' | cut -f3 > ${file/.vcf/.vcf.MAF0.01.recode.vcf.SNPs}

# Make random SNP-lists to generate differently sized samples (different numbers of SNPs included)
cat ${file/.vcf/.vcf.MAF0.01.recode.vcf.SNPs} | shuf -n100 > 100SNPs.list
cat ${file/.vcf/.vcf.MAF0.01.recode.vcf.SNPs} | shuf -n1000 > 1kSNPs.list
cat ${file/.vcf/.vcf.MAF0.01.recode.vcf.SNPs} | shuf -n10000 > 10kSNPs.list
cat ${file/.vcf/.vcf.MAF0.01.recode.vcf.SNPs} | shuf -n100000 > 100kSNPs.list

# Extract those samples
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --snps 100SNPs.list --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.100SNPs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --snps 1kSNPs.list --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.1kSNPs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --snps 10kSNPs.list --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.10kSNPs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --snps 100kSNPs.list --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.100kSNPs} &>/dev/null

# Remove the SNP-lists
rm *.list

# Extract samples of different sizes (different numbers of individuals included)
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --max-indv 25 --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.25INDs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --max-indv 100 --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.100INDs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --max-indv 500 --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.500INDs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --max-indv 1000 --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.1000INDs} &>/dev/null
vcftools --vcf ${file/.vcf/.vcf.MAF0.01.recode.vcf} --max-indv 2500 --recode --out ${file/.vcf/.vcf.MAF0.01.recode.vcf.2500INDs} &>/dev/null

# Generate IMPUTE files for each sample
for vcf in *s.recode.vcf; do vcftools --vcf $vcf --IMPUTE --out $vcf &>/dev/null; done

# Remove unnecessary files
rm *.log
rm *.legend
rm *.indv

# Generate a genetic map-file for each sample
for vcf in *s.recode.vcf; do cat $vcf | grep -v '#' | awk '{ print $1 "\t" $3 "\t" $2 "\t"  $2}' > ${vcf/s.recode.vcf/s.recode.vcf.rmap}; done

# Benchmark hapbin
for vcf in *s.recode.vcf; do taskset 1 /usr/bin/time -v ../tools/hapbin/build/ihsbin --hap ${vcf/s.recode.vcf/s.recode.vcf.impute.hap} --map ${vcf/s.recode.vcf/s.recode.vcf.rmap} --minmaf 0.0 --out '../results/benchmark/'${vcf/s.recode.vcf/s.recode.vcf.impute.hap.out} &> '../results/benchmark/'${vcf/s.recode.vcf/s.recode.vcf.impute.hap.out.bench}; done

# Extract results for individuals
for file in ../results/benchmark/*INDs.*.out.bench; do (cat $file | grep 'Haplotype count' | cut -d' ' -f3; cat $file | grep 'User time' | cut -d' ' -f4; cat $file | grep 'Maximum resident' | cut -d' ' -f6) > ${file/.bench/.benchN}; done
paste ../results/benchmark/*INDs.*.out.benchN | awk '{ print $4 "\t" $2 "\t" $5 "\t" $1 "\t" $3}' > ../results/benchmark/hapbinIND.bench

# Extract results for SNPs
for file in ../results/benchmark/*SNPs.*.out.bench; do (cat $file | grep 'Loaded' | cut -d' ' -f2; cat $file | grep 'User time' | cut -d' ' -f4; cat $file | grep 'Maximum resident' | cut -d' ' -f6) > ${file/.bench/.benchN}; done
paste ../results/benchmark/*SNPs.*.out.benchN | awk '{ print $2 "\t" $4 "\t" $3 "\t" $1}' > ../results/benchmark/hapbinSNP.bench

# Remove temporary files
rm ../results/benchmark/*.out.benchN

# Compute lookup-table for 1000 Genomes phase 3
../tools/haploblocks/filter_lookup -max_k 5008 > 1000GPP3.lookup

# Benchmark haploblocks
for vcf in *s.recode.vcf; do taskset 1 /usr/bin/time -v ../tools/haploblocks/full --out_folder ./ --vcf_path $vcf --genetic_map_path plink.chr22.GRCh37.map --lookup_path 1000GPP3.lookup --remove &> '../results/benchmark/'${vcf/s.recode.vcf/s.recode.vcf.blocks.bench}; done

# Extract results for individuals
for file in ../results/benchmark/*INDs.*.blocks.bench; do (cat $file | grep 'individuals' | cut -d' ' -f1; cat $file | grep 'User time' | cut -d' ' -f4; cat $file | grep 'Maximum resident' | cut -d' ' -f6) > ${file/.bench/.benchN}; done
paste ../results/benchmark/*INDs.*.blocks.benchN | awk '{ print $4 "\t" $2 "\t" $5 "\t" $1 "\t" $3}' > ../results/benchmark/haploblocksIND.bench

# Extract results for SNPs
for file in ../results/benchmark/*SNPs.*.blocks.bench; do (cat $file | grep 'Wrote' | cut -d' ' -f2; cat $file | grep 'User time' | cut -d' ' -f4; cat $file | grep 'Maximum resident' | cut -d' ' -f6) > ${file/.bench/.benchN}; done
paste ../results/benchmark/*SNPs.*.blocks.benchN | awk '{ print $2 "\t" $4 "\t" $3 "\t" $1}' > ../results/benchmark/haploblocksSNP.bench

# Remove temporary files
rm ../results/benchmark/*.blocks.benchN

# Plot results
Rscript ../scripts/Plot_benchmark.R ../results/benchmark/

cd ../
