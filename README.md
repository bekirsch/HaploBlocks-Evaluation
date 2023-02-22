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

## Figure 1 a)

<details>
    <summary>Show code</summary>

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

   3. Run 50 simulations (make use of GNU parallel to speed up - provided your setup allows it), one of which is started via:
   ```
   for i in {1..50}; do simulating; done
   ```

   4. Create a lookup-table:
   ```
   tools/haploblocks/filter_lookup -max_k 2000 > results/evaluation/Additive_10Mb_10kNe/ancestry.lookup
   ```

   5. Create a directory for the output:
   ```
   mkdir results/evaluation/Additive_10Mb_10kNe/output
   ```

   6. Define a function for running haploblocks:
   ```
   haploblocks() {
   vcf_gz=$1;
   vcf=${vcf_gz/.vcf.gz/.vcf}
   cmap=${vcf/.vcf/.vcf.positions};
   rmap=${cmap/.positions/.positions.map};

   zcat $vcf_gz > $vcf
   tools/haploblocks/extract_positions -i $vcf -o $cmap &>/dev/null;
   awk -v OFS='\t' '{print "chr1", "snp"NR, (50*log(1/(1-(2*1e-8*$0)))), $0}' $cmap | tr ',' '.' > $rmap;
   tools/haploblocks/full --out_folder results/evaluation/Additive_10Mb_10kNe/output --vcf_path $vcf --genetic_map_path $rmap --lookup_path results/evaluation/Additive_10Mb_10kNe/ancestry.lookup --remove &>/dev/null;

   rm $vcf;
   rm $cmap;
   rm $rmap;
   }
   export -f haploblocks
   ```

   7. Run haploblocks:
   ```
   for file in results/evaluation/Additive_10Mb_10kNe/simulation*/*.uniform.vcf.gz; do haploblocks $file; done
   ```
   8. Count the simulations (needed for plotting):
   ```
   a
   9. Plot Figure:
   ```
   Rscript scripts/Plot_Fig1a.R results/evaluation/Additive_10Mb_10kNe/output $files
   ```ll=$(ls results/evaluation/Additive_10Mb_10kNe/output/*filtered.sHat.csv | wc -l)
   files=$((all / 13))
   ```

   9. Plot Figure:
   ```
   Rscript scripts/Plot_Fig1a.R results/evaluation/Additive_10Mb_10kNe/output $files
   ```

</details>

## Figure 1 b)


### Software Versions:

SLiM:       3.4

tskit:      0.3.2

pyslim:     0.403

msprime:    0.7.4

<!---

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

3. Run 50 simulations (make use of GNU parallel to speed up - provided your setup allows it), one of which is started via:
```
for i in {1..50}; do simulating; done
```

4. Create a lookup-table:
tools/haploblocks/filter_lookup -max_k 2000 > results/evaluation/```
```
Additive_10Mb_10kNe/ancestry.lookup
```

5. Create a directory for the output:
```
mkdir results/evaluation/Additive_10Mb_10kNe/output
```

6. Define a function for running haploblocks:
```
haploblocks() {
vcf_gz=$1;
vcf=${vcf_gz/.vcf.gz/.vcf}
cmap=${vcf/.vcf/.vcf.positions};
rmap=${cmap/.positions/.positions.map};

zcat $vcf_gz > $vcf
tools/haploblocks/extract_positions -i $vcf -o $cmap &>/dev/null;
awk -v OFS='\t' '{print "chr1", "snp"NR, (50*log(1/(1-(2*1e-8*$0)))), $0}' $cmap | tr ',' '.' > $rmap;
tools/haploblocks/full --out_folder results/evaluation/Additive_10Mb_10kNe/output --vcf_path $vcf --genetic_map_path $rmap --lookup_path results/evaluation/Additive_10Mb_10kNe/ancestry.lookup --remove &>/dev/null;

rm $vcf;
rm $cmap;
rm $rmap;
}
export -f haploblocks
```

7. Run haploblocks:
```
for file in results/evaluation/Additive_10Mb_10kNe/simulation*/*.uniform.vcf.gz; do haploblocks $file; done
```
8. Count the simulations (needed for plotting):
```
a
9. Plot Figure:
```
Rscript scripts/Plot_Fig1a.R results/evaluation/Additive_10Mb_10kNe/output $files
```ll=$(ls results/evaluation/Additive_10Mb_10kNe/output/*filtered.sHat.csv | wc -l)
files=$((all / 13))
```

9. Plot Figure:
```
Rscript scripts/Plot_Fig1a.R results/evaluation/Additive_10Mb_10kNe/output $files
```



# Run selection scan

