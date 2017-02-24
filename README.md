# About
This repository contains code supporting the paper _Coordinates and Intervals on Offset-based graphs_.
See the section "How to run the gene experiments" for information about how to run the gene experiments.

## Web-tool for visualising genes on GRCh38

For a demo of the interactive web tool, please follow this link:  http://46.101.93.163/gen-graph-coords/

## Requirements
There are two requirements:
* The python package _offsetbasedgraph_. Install it with `pip3 install offsetbasedgraph`.
* A small python-package for fetching genomic sequences from UCSC. Install it by cloning the repo and installing with pip:
```
git clone git@github.com:uio-cels/gendatafetcher.git
cd gendatafetcher
pip3 install -e .
```

# How to run the gene experiments

## Analysing the


This repository contains a collection of scripts that uses the python package _offsetbasedgraph_ to create sample graphs from GRCh38.
It does so by merging flanking regions of the alt loci with the main path.

# Run example

 To run from commandline:
 
> python interface.py align_region_html chrN_GENBANKSEQACCS_alt

Example:

> python interface.py align_region_html chr7_KI270808v1_alt

Note: These scripts return html that is interpreted in web-gui/index.html.
The current demo runs web-gui/index.html and uses the php-script web-gui/python_runner.py to
call the python script Interface.py.



...