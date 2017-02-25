"""
Interace for interacting with the package OffsetBasedGraphs

Usage example:
$ python3 gen_graph_coords.py grch38.chrom.sizes grch38_alt_loci.txt genes_chr1_GL383518v1_alt.txt

With only two alt loci:
python3 gen_graph_coords.py grch38.chrom.sizes-small grch38_alt_loci_small.txt genes_chr1_GL383518v1_alt.txt

"""

from collections import defaultdict
import sys
import argparse


from methods import *

# Dict struct for holding all arguments taken by the interface

CHROM_SIZES_DESCRIPTION = 'Name of file with chrom sizes, used to build graph (e.g. data/grch38.chrom.sizes).' + \
                          ' Should contain two columns, chrom/alt name and size.'


interface = \
{
    'create_graph':
        {
            'help': 'Create graph',
            'arguments':
                [
                    ('chrom_sizes_file_name', CHROM_SIZES_DESCRIPTION),
                    ('alt_locations_file_name', 'File containing alternative loci info (e.g. data/grch38_alt_loci.txt)'),
                    ('out_file_name', 'Name of file to store graph and translation objects insize')
                ],
            'method': create_graph
        },

    'check_duplicate_genes':
        {
            'help': 'Experiment: Analyse duplicate genes on graph create from GRCh38',
            'arguments':
                [
                    ('translation_file_name', 'Translation file created by running create_graph'),
                    ('genes_file_name', '')
                ],
            'method': check_duplicate_genes
        },

    'create_complex_graph':
        {
            'help': 'Create a complex graph from GRCh38 by merging alt loci in using NCBI alignments.',
            'arguments':
                [
                    ('chrom_sizes_file_name', CHROM_SIZES_DESCRIPTION),
                    ('out_file_name', 'File to store resulting translation object in')
                ],
            'method': merge_all_alignments
        },

    'analyse_multipath_genes':
        {
            'help': 'Experiment: Analyse multi-path genes on a graph created from GRCh38\n' +
                    'Run example: \n $ python3 gen_graph_coords.py analyse_multipath_genes' +
                    'data/grch38.chrom.sizes data/grch38_alt_loci.txt data/alt_alignments' +
                    'data/genes/genes_refseq.txt',
            'arguments':
                [
                    ('chrom_sizes_file_name', CHROM_SIZES_DESCRIPTION),
                    ('alt_locations_file_name', 'File containing alternative loci info (e.g. data/grch38_alt_loci.txt)'),
                    ('ncbi_alignments_dir', 'Directory containing NCBI alignment files (e.g. data/alt_alignments)'),
                    ('genes_file_name', 'Name of gene file (e.g. data/genes/genes_refseq.txt)'),
                ],
            'method': analyse_multipath_genes2
        },
    'visualize_alt_locus':
        {
            'help': 'Produce html visualization (that can be saved and opened in a browser)',
            'arguments':
                [
                    ('translation_file_name', 'Translation file created by running create_graph. Contains the graph that will be visualized'),
                    ('genes', 'Name of gene file containing genes that will be visualized (e.g. data/genes/genes_refseq.txt)'),
                    ('alt_locations_file_name', 'File containing alternative loci info (e.g. data/grch38_alt_loci.txt)'),
                    ('alt_locus', 'Alt locus id (e. g. chr2_KI270774v1_alt)')
                ],
            'method': visualize_alt_locus
        },

    'visualize_alt_locus_wrapper':
        {
            'help': 'Wrapper for visualize_alt_locus. Requires specific files in specific places. Not recomended to run.',
            'arguments':
                [
                    ('alt_locus', 'Alt locus id (e. g. chr2_KI270774v1_alt)')
                ],
            'method': visualize_alt_locus_wrapper
        },

    'html_alt_loci_select':
        {
            'help': 'Produce html for alt loci select box (only used by web tool)',
            'arguments':
                [],
            'method': html_alt_loci_select
        }


}

# Create parser
parser = argparse.ArgumentParser(
    description='Interact with a graph created from GRCh38')
subparsers = parser.add_subparsers(help='Subcommands')

for command in interface:
    subparser = subparsers.add_parser(command, help=interface[command]["help"])
    for argument, help in interface[command]["arguments"]:
        subparser.add_argument(argument, help=help)
    subparser.set_defaults(func=interface[command]["method"])

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
if hasattr(args, 'func'):
    args.func(args)
else:
    parser.help()




sys.exit()