from offsetbasedgraph import Graph, Translation
from offsetbasedgraph.gene import GeneList
import sys

from offsetbasedgraph.graphutils import *


def create_graph(args):
    graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    new_numeric_graph, numeric_translation = connect_without_flanks(
        numeric_graph, args.alt_locations_file_name, name_translation)
    name_graph, new_name_translation = convert_to_text_graph(
        new_numeric_graph, name_translation, numeric_translation)
    final_translation = name_translation + numeric_translation + new_name_translation
    final_translation.graph2 = name_graph

    final_translation.to_file(args.out_file_name)
    print("Graph and translation object stored in %s" % (args.out_file_name))


def check_duplicate_genes(args):
    genes_file_name = args.genes_file_name
    final_trans = Translation.from_file(args.translation_file_name)
    genes = get_gene_objects_as_intervals(genes_file_name, final_trans.graph1)
    analyze_genes_on_merged_graph(genes, final_trans)
    # print(genes_file_name)


def merge_alignment(args):
    # For every alt loci, create complex graph,translate genes and analyse them
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    graph, name_trans = grch38_graph_to_numeric(text_graph)
    alt_locus = args.alt_locus_id
    trans, complex_graph = merge_alt_using_cigar(graph, name_trans, alt_locus)

    full_trans = name_trans + trans
    full_trans.to_file(args.out_file_name)

    # Read genes and translate to graph
    genes = get_gene_objects_as_intervals(args.genes)
    genes_on_alt = []
    for g in genes:
        if g.transcription_region.start_position.region_path_id == alt_locus:
            genes_on_alt.append(g.translate(full_trans))

    genes_on_alt = GeneList(genes_on_alt)
    genes_on_alt.to_file("genes_%s" % args.out_file_name)

    full_trans.to_file(args.out_file_name)
    print("Saved trans to file %s" % args.out_file_name)


def merge_all_alignments(args):
    from offsetbasedgraph.graphutils import merge_alt_using_cigar, \
        grch38_graph_to_numeric
    # Text ids (chrom names and alt names)
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    graph, name_trans = grch38_graph_to_numeric(text_graph)

    # Go through all alts in this graph
    new_graph = graph.copy()
    i = 0
    for b in text_graph.blocks:
        if "alt" in b:
            print("Merging %s" % b)
            numeric_trans, new_graph = merge_alt_using_cigar(
                new_graph, name_trans, b)
            i += 1
            if i >= 1:
                break

    final_graph, trans2 = convert_cigar_graph_to_text(
        new_graph,
        name_trans,
        numeric_trans)

    trans_a = name_trans + numeric_trans
    full_trans = trans_a + trans2

    # full_trans = name_trans + numeric_trans + trans2
    full_trans.set_graph2(final_graph)
    print("To file")
    full_trans.to_file(args.out_file_name)


def visualize_alt_locus_wrapper(args, quiet=False):

    if not quiet:
        print("<div style='display: none'>")

    # Finds correct gene file etc
    chrom = args.alt_locus.split("_")[0]
    args.genes = "data/genes/genes_refseq_%s.txt" % (chrom)

    # Create graph only for this alt loci
    graph = create_initial_grch38_graph("data/grch38.chrom.sizes")
    #graph.to_file("graph_non_connected")
    #graph = Graph.from_file("graph_non_connected")

    #print("Convert to numeric")
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    #name_translation.to_file("name_trans_non_connected")

    #print("Convert without flanks")
    new_numeric_graph, numeric_translation = connect_without_flanks(
        numeric_graph, "data/grch38_alt_loci.txt", name_translation, [args.alt_locus])

    name_graph, new_name_translation = convert_to_text_graph(
        new_numeric_graph, name_translation, numeric_translation)

    final_translation = name_translation + numeric_translation + new_name_translation
    final_translation.graph2 = name_graph
    #final_translation.to_file("tmp_trans")

    args.translation_file_name = final_translation
    args.alt_locations_file_name = 'data/grch38_alt_loci.txt'

    if not quiet:
        print("</div>")
    #return
    visualize_alt_locus(args, True, quiet)


def visualize_alt_locus(args, skip_wrapping=False, quiet=False):
    from offsetbasedgraph.graphutils import GeneList, \
        create_gene_dicts, create_subgraph_around_alt_locus

    if not isinstance(args.translation_file_name, Translation):
        trans = Translation.from_file(args.translation_file_name)
    else:
        trans = args.translation_file_name

    graph = trans.graph2
    orig_trans = trans.copy()

    # Find all genes on this graph
    genes = GeneList(get_gene_objects_as_intervals(args.genes)).gene_list
    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(genes, alt_loci_fn=args.alt_locations_file_name)
    genes = main_genes[args.alt_locus]
    genes = [g.translate(trans) for g in genes]
    subgraph, trans, start_position = create_subgraph_around_alt_locus(graph, trans, args.alt_locus, 200000, alt_loci_fn=args.alt_locations_file_name)

    start_position = orig_trans.translate_position(start_position, True)[0]

    genes = [g for g in genes if not g.multiple_alt_loci() and g.transcription_region.length() > 100]

    if len(genes) > 40:
        genes = genes[0:40]

    levels = Graph.level_dict(subgraph.blocks)

    # Find start block by choosing a block having no edges in
    start = None
    for b in subgraph.blocks:
        if len(subgraph.reverse_adj_list[b]) == 0:
            start = b
            break

    assert start is not None

    from visualizehtml import VisualizeHtml
    subgraph.start_block = start
    max_offset = sum([subgraph.blocks[b].length() for b in subgraph.blocks])
    v = VisualizeHtml(subgraph, 0, max_offset, 0, levels, "", 800, genes, start_position)

    if quiet:
        return

    if skip_wrapping:
        print(str(v))
    else:
        print(v.get_wrapped_html())


def _analyse_multipath_genes_on_graph(genes_list, genes_against, graph):
    # Takes a list of mp genes and a graph
    # Returns number of equal exons and equal genes
    equal = 0
    equal_exons = 0
    n = 1
    for g in genes_list:
        if n % 1000 == 0:
            print("Checked %d genes" % n)
        n += 1

        for g2 in genes_against:
            if g is g2:
                continue

            if g == g2:
                equal += 1

            if g.faster_equal_critical_intervals(g2):
                equal_exons += 1

    return equal, equal_exons


def analyze_fuzzy_genes(args):
    genes = get_gene_objects_as_intervals(args.genes_file_name)
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    return fuzzy_gene_analysis(genes, text_graph, args.ncbi_alignments_dir,
                               args.alt_locations_file_name)


def analyse_multipath_genes2(args):
    if args.interval_type == "fuzzy":
        return analyze_fuzzy_genes(args)
    assert args.interval_type == "critical"
    from offsetbasedgraph.graphutils import create_gene_dicts, \
        translate_single_gene_to_aligned_graph
    print("Reading genes")
    genes = get_gene_objects_as_intervals(args.genes_file_name)

    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(
        genes,
        alt_loci_fn=args.alt_locations_file_name)

    # alt loci genes are only genes on alt loci (nothing on main)
    # exon_dict contains only genes on main, index by offset of first exon

    # For every alt loci, create complex graph, translate genes and analyse them
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    graph, name_trans = grch38_graph_to_numeric(text_graph)

    equal_total = 0
    equal_exons_total = 0
    n_a = 1
    print(text_graph.blocks.keys())
    for b in text_graph.blocks:
        if "alt" in b:
            sys.stdout.flush()

            n_a += 1
            genes_here = alt_loci_genes[b]
            if not (genes_here and main_genes[b]):
                print("Skipping", b)
                continue

            trans, complex_graph = merge_alt_using_cigar(graph, name_trans, b, ncbi_alignments_dir=args.ncbi_alignments_dir)
            full_trans = name_trans + trans

            # Find candidates on main path to check against:
            genes_against = [g.copy() for g in main_genes[b]]
            genes_against_translated = []
            n = 0
            for mg in genes_against:
                n += 1
                genes_against_translated.append(translate_single_gene_to_aligned_graph(mg, full_trans).interval)
                sys.stdout.write('\r  Translating main genes: ' + str(round(100 * n / len(genes_against))) + ' % finished ' + ' ' * 20)
                sys.stdout.flush()

            print()

            genes_here_translated = []
            n = 0
            for mg in genes_here:
                n += 1
                genes_here_translated.append(translate_single_gene_to_aligned_graph(mg, full_trans).interval)
                sys.stdout.write('\r  Translating alt genes: ' + str(round(100 * n / len(genes_here))) + ' % finished ' + ' ' * 20)
                sys.stdout.flush()
            print()
            equal, equal_exons = _analyse_multipath_genes_on_graph(
                genes_here_translated,
                genes_against_translated,
                complex_graph)
            equal_total += equal
            equal_exons_total += equal_exons

    print("SUM:")
    print("Equal: %d, equal exons: %d" % (equal_total, equal_exons_total))

    print("RESULTS:")
    print(" Number of genes originating from alt loci identical to a gene originaing from main chromosome: %d" % equal_total)
    print(" Number of genes with only identical exones (not start and end position): %d" % equal_exons_total)


def html_alt_loci_select(args):
    # Prints all regions as an html select field
    html_out = """<select name='region'
               class='form-control' style='width: 320px;'>"""
    from offsetbasedgraph.graphutils import get_alt_loci_positions
    loci = get_alt_loci_positions("data/grch38_alt_loci.txt")

    for locus in loci:
        region = loci[locus]
        if region["end"] - region["start"] < 4000000:
            html_out += "<option value='%s'>%s (%s:%d-%d)</option>" % \
                        (locus,
                         locus,
                         region["main_chr"], \
                         region["start"], region["end"])
    html_out += "</select>"
    print(html_out)


def print_gene_notations(args):
    print("Processing genes")
    genes_file_name = args.genes
    trans = Translation.from_file(args.translation_file_name)


    genes = get_gene_objects_as_intervals(genes_file_name, trans.graph1)
    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(genes, alt_loci_fn=args.alt_locations_file_name)
    genes = main_genes[args.alt_locus]

    print("Original genes on GRCh38")
    print("----------------------------")
    for g in genes:
        print("%s: %s" % (g.name, g.transcription_region.notation()))

    print("Genes on graph")

    genes = [g.translate(trans) for g in genes]

def compute_average_flank_length(args):
    from offsetbasedgraph.GRCH38 import AltLoci
    import numpy as np

    alts = AltLoci.from_file(args.alt_locations_file_name)
    lengths = []
    for alt in alts.alt_loci:
        lengths.append(alt.start_flank.length())
        lengths.append(alt.end_flank.length())

    print(np.mean(lengths))

