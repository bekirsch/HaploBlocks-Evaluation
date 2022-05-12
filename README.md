# Evaluating haploblocks

This repository contains all code used for evaluating [HaploBlocks](https://github.com/bekirsch/HaploBlocks) on simulated data as well as an example on how to perform a chromosome-wide scan.

We assume from here on that you have cloned this github repository into a directory called e.g.
`haploblocks-evaluation`, and are running commands from within it, e.g. using

```
git clone https://github.com/BeGerweck/haploblocks-evaluation.git
cd haploblocks-evaluation
```

# Requirements and dependencies

All of the below commands were tested on Ubuntu 20.04 LTS. If you haven't already, download and make haploblocks from https://github.com/bekirsch/HaploBlocks, preferably in `haploblocks-evaluation/tools`. 

```
cd tools
git clone https://github.com/bekirsch/HaploBlocks
cd haploblocks
make
cd ../../
```

## Install prerequisites

You will need Python (version 3) with pip and the GNU scientific library (`gsl`) (for msprime/pyslim), as well as R (for plotting) and cmake (for SLiM). To install these on Ubuntu:

```
sudo apt-get install python3-pip libgsl-dev r-base-core cmake
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

Download and install [SLiM](http://messerlab.org/slim/) following their user manual. Note, that the following simulations were conducted in SLiM 3.

# Run evaluation

Make sure you have at least 40GB of available space before running the evaluation. The script `run_evaluation.sh` can be run to produce figures ... and the data needed for them. Note that this may take days or weeks depending on your system. You can run this script multithreaded, resulting in it running only hours or days, by providing a number `n` of threads. If you choose more threads than available cores, the script will run with `n = (available cores)-1`.

```
chmod +x run_evaluation.sh
./run_evaluation.sh n
```

Figures ... and the underlying data can be produced by running `script2`. This also may take days or weeks depending on your system.

```
chmod +x script2
./script2
```

ATTENTION: You can run script2 multithreaded in the same manner as run_evaluation.sh, but note that one thread in script2 can take up to 12GB of memory!

# Run selection scan
TBA
