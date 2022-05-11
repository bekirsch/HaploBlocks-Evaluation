# Evaluating haploblocks

This repository contains all code used for evaluating [haploblocks](https://gitlab.com/bacazaux/haploblocks)
and producing the figures in the paper:
[Link to paper](https://www.google.com). This
includes code for evaluating haploblocks with simulated data, benchmarking against [hapbin](https://github.com/evotools/hapbin) using real data, as well as performing a chromosome-wide scan.

We assume from here on that you have cloned this github repository into a directory called e.g.
`haploblocks-evaluation`, and are running commands from within it, e.g. using

```
git clone https://github.com/BeGerweck/haploblocks-evaluation.git
cd haploblocks-evaluation
```

# Requirements and dependencies

All of the below commands were tested on Ubuntu 20.04 LTS. If you haven't already, download and make haploblocks from https://gitlab.com/bacazaux/haploblocks, preferably in `haploblocks-evaluation/tools`. (FOR NOW THIS IS NOT NEEDED AS HAPLOBLOCKS IS ALREADY INCLUDED IN `tools`, BUT YOU MUST STILL "make" IN ORDER TO RUN HAPLOBLOCKS!)

```
cd tools
git clone https://gitlab.com/bacazaux/haploblocks
cd haploblocks
make
cd ../../
```

## Install prerequisites

You will need Python (version 3) with pip and the GNU scientific library (`gsl`) (for msprime/pyslim), as well as R (for plotting), GNU parallel (for multithreading) and cmake (for SLiM and hapbin). To install these on Ubuntu:

```
sudo apt-get install python3-pip libgsl-dev r-base-core cmake parallel
```

## Install python modules

Required Python modules are listed in `modules.txt`. Simply install them via

```
python3 -m pip install -r modules.txt
```

## Install R packages

We require the `latex2exp` and `stringr` packages. If you don't have these installed in your local R installation (make sure you have one on your system), you should be able to install them from within R via `install.packages(c("latex2exp", "stringr"))`. If you have root access to your machine, you can install packages without requiring any user interaction by
```
sudo R -e 'install.packages(c("latex2exp", "stringr"), repos="https://cran.r-project.org")'
```

## Install SLiM

Download and install [SLiM 3](http://messerlab.org/slim/) following their user manual.

## Install hapbin

### hapbin dependencies

We benchmarked our tool against [hapbin](https://github.com/evotools/hapbin), for which `vcftools` and additional dependencies for multitheading need to be installed. vcftools is needed to convert VCF files to `IMPUTE` file formats as hapbin does not accept VCF files as input.

```
sudo apt-get install vcftools mpich libmpich-dev
```
### Build hapbin

Once these dependencies are installed you can build hapbin from source code:

```
cd tools/
git clone https://github.com/evotools/hapbin.git
cd hapbin/build/
cmake ../src/
make
cd ../../../
```

# Run evaluation

Make sure you have at least 40GB of available space before running the evaluation. The script `run_evaluation.sh` can be run to produce figures ... and the data needed for them. Note that this may take days or weeks depending on your system. You can run this script multithreaded, resulting in it running only hours or days, by providing a number `n` of threads. If you choose more threads than available cores, the script will run with `n = (available cores)-1`.

```
chmod +x run_evaluation.sh
./run_evaluation.sh n
```

Figures ... and the underlying data can be produced by running `script2`. This may take days or weeks depending on your system.

```
chmod +x script2
./script2
```

ATTENTION: You can run script2 multithreaded in the same manner as run_evaluation.sh, but note that one thread in script2 can take up to 12GB of memory!

# Run benchmark

The script `run_benchmark.sh` will not produce figure??? for chromosome??? of the UK Biobank data set, but for chromsome 22 of the 1000 Genomes Project Phase 3 data set, obtainable from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, as this is pubicly available. Additionally you'll need to download genetic maps from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip. Place the gzipped VCF for chromosome 22 and the matching map `plink.chr22.GRCh37.map` in the folder `chr` before executing the script.

```
chmod +x run_benchmark.sh
./run_benchmark.sh
```

# Run selection scan
