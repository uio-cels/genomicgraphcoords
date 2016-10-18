# Live demo

For a live demo of this tool, please follow this link:  http://46.101.93.163/gen-graph-coords/

# About
This tool uses and requires the python package _offsetbasedgraph_. Install it with `pip3 install offsetbasedgraph`.

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