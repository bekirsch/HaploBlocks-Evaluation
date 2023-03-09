# Evaluating HaploBlocks

This repository contains all code used for evaluating [HaploBlocks](https://github.com/bekirsch/HaploBlocks) on simulated data as well as an example on how to perform a chromosome-wide scan.

We assume from here on that you have cloned this github repository into a directory called e.g.
`cd HaploBlocks-Evaluation`, and are running commands from within it, e.g. using

```
git clone https://github.com/bekirsch/HaploBlocks-Evaluation.git
cd HaploBlocks-Evaluation
```

# Requirements and dependencies

All of the below commands were tested on Ubuntu 20.04. If you haven't already, download and make HaploBlocks from https://github.com/bekirsch/HaploBlocks, preferably in `HaploBlocks-Evaluation/tools`. 

```
cd tools
git clone https://github.com/bekirsch/HaploBlocks
cd HaploBlocks/build
make
cd ../../../
```

## Install prerequisites

You will need Python (version 3) with pip and the GNU scientific library (`gsl`) (for msprime/pyslim), as well as R (for plotting) and cmake (for SLiM). To install these on Ubuntu:

```
sudo apt install python3-pip libgsl-dev r-base-core cmake
```

**Note:** Different simulations were performed on different versions of the dependencies below, as the programs progressed along our analyses. Specific versions will be mentioned at the relevant section of this README. Keep in mind before installing!

## Install python modules

Required Python modules are listed in `modules.txt`. You can install them via

```
python3 -m pip install -r modules.txt
```

## Install R packages

We require the `latex2exp` and `stringr` packages. If you don't have these installed in your local R installation (make sure you have one on your system), you should be able to install them from within R via `install.packages(c("latex2exp", "stringr"))`.

## Install SLiM

Download and install [SLiM](http://messerlab.org/slim/) following their user manual.

# Validation on simulated data

<details>
    <summary>Figure 1 a)</summary>

   ### Software Versions:

   SLiM:       3.4

   tskit:      0.3.2

   pyslim:     0.403

   msprime:    0.7.4

   1. Create a directory for simulations:
   ```
   mkdir results/evaluation/Additive_10Mb_10kNe
   ```
   2. Define a function for simulating:
   ```
   simulating() {
   seed=$(openssl rand 4 | od -DAn);
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
   ```

   3. Run 50 simulations (make use of GNU parallel to speed up - provided your setup allows it):
   ```
   for i in {1..50}; do simulating; done
   ```

   4. Create a lookup-table:
   ```
   tools/HaploBlocks/build/filter_lookup -max_k 2000 > results/evaluation/Additive_10Mb_10kNe/ancestry.lookup
   ```

   5. Create a directory for the output:
   ```
   mkdir results/evaluation/Additive_10Mb_10kNe/output
   ```

   6. Define a function for running HaploBlocks:
   ```
   HaploBlocks() {
   vcf_gz=$1;
   vcf=${vcf_gz/.vcf.gz/.vcf}
   cmap=${vcf/.vcf/.vcf.positions};
   rmap=${cmap/.positions/.positions.map};

   zcat $vcf_gz > $vcf
   tools/HaploBlocks/build/extract_positions -i $vcf -o $cmap &>/dev/null;
   awk -v OFS='\t' '{print "chr1", "snp"NR, (50*log(1/(1-(2*1e-8*$0)))), $0}' $cmap | tr ',' '.' > $rmap;
   tools/HaploBlocks/build/full --out_folder results/evaluation/Additive_10Mb_10kNe/output --vcf_path $vcf --genetic_map_path $rmap --lookup_path results/evaluation/Additive_10Mb_10kNe/ancestry.lookup --remove &>/dev/null;

   rm $vcf;
   rm $cmap;
   rm $rmap;
   }
   export -f HaploBlocks
   ```

   7. Run HaploBlocks:
   ```
   for file in results/evaluation/Additive_10Mb_10kNe/simulation*/*.uniform.vcf.gz; do HaploBlocks $file; done
   ```
   8. Count the simulations (needed for plotting):
   ```
   all=$(ls results/evaluation/Additive_10Mb_10kNe/output/*filtered.sHat.csv | wc -l)
   files=$((all / 13))
   ```

   9. Plot Figure:
   ```
   Rscript scripts/Plot_Fig1a.R results/evaluation/Additive_10Mb_10kNe/output $files
   ```

</details>

<details>
    <summary>Figure 1 b)</summary>

   ### Software Versions:

   SLiM:       3.4

   tskit:      0.3.2

   pyslim:     0.403

   msprime:    0.7.4

1. create directory for simulations:
```
mkdir results/evaluation/Gravel_CEU
```

2. Define a function for simulating:
```
simulating_Gravel_CEU() {
  seed=$(openssl rand 4 | od -DAn);
  slim -s $seed -d gen=$1 scripts/Gravel_CEU.slim;
}
export -f simulating_Gravel_CEU
```

3. Define a function for recapitation:
```
recap_Gravel_CEU() {
  trees=$1
  python3 scripts/recapitation_gravel_CEU.py -i $trees;
  line=$(cat ${trees/.trees/.trees.vcf} | grep -n '5000000' | cut -f1 | cut -d":" -f1)
  cat ${trees/.trees/.trees.vcf} | awk -F '\t' -v OFS='\t' -v m=$line -v n=4 -v el='0' 'NR == m { $n = el } 1' | awk -F '\t' -v OFS='\t' -v m=$line -v n=5 -v el='1' 'NR == m { $n = el } 1' | gzip > ${trees/.trees/.trees.uniform.vcf.gz}
  rm ${trees/.trees/.trees.vcf}
  line=$(zcat ${trees/.trees/.trees.uniform.vcf.gz} | cut -f2 | grep -n '\<5000000\>' | cut -d':' -f1);
  zcat ${trees/.trees/.trees.uniform.vcf.gz} | awk -v li=$line 'NR==li' | cut -f1-9 --complement | grep -o '1' | wc -l > ${trees/.trees/.trees.uniform.vcf.gz.count}
}
export -f recap_Gravel_CEU
```
4. Run simulations:

We run more simulations for intermediate generations ago, to ensure sufficient intermediate frequencies are reached.

  For upper and lower frequencies:
  ```
  for g in 5850 5800 5750 5700 5650 5300 5250 5200 5150; do
	  for i in {1..20}; do simulating_Gravel_CEU $g; done
  done
  ```
  For intermediate frequencies:
  ```
  for g in 5600 5550 5500 5450 5400 5350; do
	  for i in {1..40}; do simulating_Gravel_CEU $g; done
  done
  ```
  And to recapitate:
  ```
  for file in results/evaluation/Gravel_CEU/*simulation*/*trees; do recap_Gravel_CEU $file; done
  ```
5. Create directory for output:
```
mkdir results/evaluation/output_gravel_CEU
mv results/evaluation/Gravel_CEU/*simulation*/*.count results/evaluation/output_gravel_CEU
```
6. Define a function for running HaploBlocks:
```
HaploBlocks_gravel_CEU() {
vcf_gz=$1;
vcf=${vcf_gz/.vcf.gz/.vcf}
cmap=${vcf/.vcf/.vcf.positions};
rmap=${cmap/.positions/.positions.map};
zcat $vcf_gz > $vcf
tools/HaploBlocks/extract_positions -i $vcf -o $cmap &>/dev/null;
awk -v OFS='\t' '{print "chr1", "snp"NR, (50*log(1/(1-(2*1e-8*$0)))), $0}' $cmap | tr ',' '.' > $rmap;
tools/HaploBlocks/full --out_folder results/evaluation/output_gravel_CEU --vcf_path $vcf --genetic_map_path $rmap --lookup_path results/evaluation/ancestry.lookup --remove &>/dev/null;
rm $vcf;
rm $cmap;
rm $rmap;
}
export -f HaploBlocks_gravel_CEU
```
7. Run HaploBlocks:
```
for file in results/evaluation/Gravel_CEU/*simulation*/*.uniform.vcf.gz; do HaploBlocks_gravel_CEU $file; done
```

8. Plot results:
```
Rscript scripts/Plot_Fig1b.R results/evaluation/output_gravel_CEU $files
```


</details>

# Benchmarks

<details>
    <summary>Figure 2</summary>

   ### Software Versions:

   RAISD:       2.9


1. Download and install [RAiSD](https://github.com/alachins/raisd) to the following their instructions, again preferably in `HaploBlocks-Evaluation/tools`.

2. Download and unpack [chromosome 2](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz):

**Note:** The unpacked size of this vcf is ~72GB. You can use a smaller chromosome instead.

```
mkdir chr
cd chr/
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
gzip -d ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

3. Download the appropriate genetic maps, i.e.:
```
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
unzip plink.GRCh37.map.zip
```

4. Install vcftools:
```
sudo apt install vcftools
```

5. Remove indels and for MAF > 0.01:
```
file=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf
vcftools --vcf $file --remove-indels --maf 0.01 --recode --out ${file/.vcf/.MAF0.01}
```

6. Extract a list of all SNPs:
```
file=ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.MAF0.01.recode.vcf
cat $file | grep -v '#' | cut -f3 > ${file/.vcf/.SNPs}
```

7. Generate random SNP-lists:
```
cat ${file/.vcf/.SNPs} | shuf -n100 > 100SNPs.list
cat ${file/.vcf/.SNPs} | shuf -n1000 > 1kSNPs.list
cat ${file/.vcf/.SNPs} | shuf -n10000 > 10kSNPs.list
cat ${file/.vcf/.SNPs} | shuf -n100000 > 100kSNPs.list
```

8. Generate vcfs containing only these SNPs:
```
vcftools --vcf $file --snps 100SNPs.list --recode --out ${file/.recode.vcf/.100SNPs}
vcftools --vcf $file --snps 1kSNPs.list --recode --out ${file/.recode.vcf/.1kSNPs}
vcftools --vcf $file --snps 10kSNPs.list --recode --out ${file/.recode.vcf/.10kSNPs}
vcftools --vcf $file --snps 100kSNPs.list --recode --out ${file/.recode.vcf/.100kSNPs}
```

9. Generate vcfs containg different numbers of random individuals:
```
vcftools --vcf $file --max-indv 5 --recode --out ${file/.recode.vcf/.5INDs}
vcftools --vcf $file --max-indv 50 --recode --out ${file/.recode.vcf/.50INDs}
vcftools --vcf $file --max-indv 500 --recode --out ${file/.recode.vcf/.500INDs}
```

10. Run RAiSD on all vcfs:
```
for vcf in *.recode.vcf; do ../tools/RAiSD -n $vcf -I $vcf -f; done
```

</details>

<!---

 Generate a genetic map-file for each sample
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

--->