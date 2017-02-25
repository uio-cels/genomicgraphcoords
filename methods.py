from offsetbasedgraph import Graph, Translation
from offsetbasedgraph.gene import GeneList, Gene, MultiPathGene

from offsetbasedgraph.graphutils import Gene, convert_to_numeric_graph, connect_without_flanks, \
    convert_to_text_graph, merge_flanks, connect_without_flanks, parse_genes_file, \
    get_genes_as_intervals, get_gene_objects_as_intervals, find_exon_duplicates, \
    create_initial_grch38_graph, blast_test, convert_cigar_graph_to_text, analyze_genes_on_merged_graph, \
    merge_alt_using_cigar, grch38_graph_to_numeric, create_subgraph_around_alt_locus

def create_graph(args):
    graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    n_starts = len([b for b in graph.blocks if not graph.reverse_adj_list[b]])
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    n_starts2 = len([b for b in numeric_graph.blocks if not numeric_graph.reverse_adj_list[b]])
    assert n_starts == n_starts2
    new_numeric_graph, numeric_translation = connect_without_flanks(
        numeric_graph, args.alt_locations_file_name, name_translation)
    n_starts3 = len([b for b in new_numeric_graph.blocks if not new_numeric_graph.reverse_adj_list[b]])
    name_graph, new_name_translation = convert_to_text_graph(
        new_numeric_graph, name_translation, numeric_translation)
    n_starts4 = len([b for b in name_graph.blocks if not name_graph.reverse_adj_list[b]])
    assert n_starts3 == n_starts4
    final_translation = name_translation + numeric_translation + new_name_translation
    final_translation.graph2 = name_graph
    n_starts5 = len([b for b in final_translation.graph2.blocks if not final_translation.graph2.reverse_adj_list[b]])

    assert n_starts5 == n_starts4
    final_translation.to_file(args.out_file_name)
    print("Graph and translation object stored in %s" % (args.out_file_name))


def check_duplicate_genes(args):
    genes_file_name = args.genes_file_name
    final_trans = Translation.from_file(args.translation_file_name)
    genes = get_gene_objects_as_intervals(genes_file_name, final_trans.graph1)
    analyze_genes_on_merged_graph(genes, final_trans)
    # print(genes_file_name)


def merge_alignment(args):
    from offsetbasedgraph.graphutils import GeneList, create_gene_dicts, translate_single_gene_to_aligned_graph

    # For every alt loci, create complex graph, translate genes and analyse them
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    graph, name_trans = grch38_graph_to_numeric(text_graph)
    alt_locus = args.alt_locus_id
    trans, complex_graph = merge_alt_using_cigar(graph, name_trans, alt_locus)

    full_trans = name_trans + trans
    full_trans.to_file(args.out_file_name)

    # Read genes and translate to graph
    genes = GeneList(get_gene_objects_as_intervals(args.genes)).gene_list
    genes_on_alt = []
    for g in genes:
        if g.transcription_region.start_position.region_path_id == alt_locus:
            genes_on_alt.append(g.translate(full_trans))

    genes_on_alt = GeneList(genes_on_alt)
    genes_on_alt.to_file("genes_%s" % args.out_file_name)
    #print("Genes on alt")
    #print(genes_on_alt)

    full_trans.to_file(args.out_file_name)
    print("Saved trans to file %s" % args.out_file_name)

def merge_all_alignments(args):
    from offsetbasedgraph.graphutils import merge_alt_using_cigar, grch38_graph_to_numeric
    # Text ids (chrom names and alt names)
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    graph, name_trans = grch38_graph_to_numeric(text_graph)

    # Go through all alts in this graph
    new_graph = graph.copy()
    i = 0
    #for b in text_graph.blocks:
    for b in ['chr8_KI270818v1_alt']: #['chr8_KI270812v1_alt']: #text_graph.blocks: # chr6_GL000251v2_alt

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


def visualize_alt_locus_wrapper(args):
    # Finds correct gene file etc
    chrom = args.alt_locus.split("_")[0]
    args.genes = "genes/genes_refseq_%s.txt" % (chrom)

    # Create graph only for this alt loci
    graph = create_initial_grch38_graph("grch38.chrom.sizes")
    #graph.to_file("graph_non_connected")
    #graph = Graph.from_file("graph_non_connected")

    #print("Convert to numeric")
    numeric_graph, name_translation = convert_to_numeric_graph(graph)
    #name_translation.to_file("name_trans_non_connected")

    #print("Convert without flanks")
    new_numeric_graph, numeric_translation = connect_without_flanks(
        numeric_graph, "grch38_alt_loci.txt", name_translation, [args.alt_locus])

    name_graph, new_name_translation = convert_to_text_graph(
        new_numeric_graph, name_translation, numeric_translation)

    final_translation = name_translation + numeric_translation + new_name_translation
    final_translation.graph2 = name_graph
    #final_translation.to_file("tmp_trans")

    args.translation_file_name = final_translation
    #return
    visualize_alt_locus(args, True)


def visualize_alt_locus(args, skip_wrapping=False):
    from offsetbasedgraph.graphutils import GeneList, create_gene_dicts, create_subgraph_around_alt_locus

    if not isinstance(args.translation_file_name, Translation):
        trans = Translation.from_file(args.translation_file_name)
    else:
        trans = args.translation_file_name

    graph = trans.graph2
    orig_trans = trans.copy()


    # Find all genes on this graph
    genes = GeneList(get_gene_objects_as_intervals(args.genes)).gene_list

    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(genes)

    alt = args.alt_locus
    genes = alt_loci_genes[args.alt_locus] + main_genes[args.alt_locus]
    genes = main_genes[args.alt_locus]

    #print("Number of genes: %d" % (len(genes)))
    genes = [g.translate(trans) for g in genes]
    trans_regions = [g.transcription_region for g in genes]

    if len(trans_regions) == 0:
        raise Exception("No genes in area")

    #subgraph, trans, start_position = graph.create_subgraph_from_intervals(trans_regions, 200000, args.alt_locus)
    subgraph, trans, start_position = create_subgraph_around_alt_locus(graph, trans, args.alt_locus, 200000)

    start_position = orig_trans.translate_position(start_position, True)[0]
    print("<p>Start position: %s</p>" % start_position)
    #trans.graph1 = full_trans.graph2

    #full_trans = full_trans + trans

    genes = [g.translate(trans) for g in genes]

    genes = [g for g in genes if not g.multiple_alt_loci()]


    if len(genes) > 3:
        #genes.sort(key=lambda g: g.length(), reverse=True)
        genes = genes[0:3]

    levels = Graph.level_dict(subgraph.blocks)

    # Find start block by choosing a block having no edges in
    start = None
    for b in subgraph.blocks:
        if len(subgraph.reverse_adj_list[b]) == 0:
            start = b
            break

    #print("== Subgraph ==")
    #print(subgraph)

    assert start is not None

    #full_trans = full_trans + trans

    from offsetbasedgraph import VisualizeHtml
    subgraph.start_block = start
    max_offset = sum([subgraph.blocks[b].length() for b in subgraph.blocks])
    v = VisualizeHtml(subgraph, 0, max_offset, 0, levels, "", 800, genes, start_position)

    if skip_wrapping:
        print(str(v))
    else:
        print(v.get_wrapped_html())


def visualize_genes(args):
    trans = Translation.from_file(args.translation_file_name)
    graph = trans.graph2

    # Find blocks that genes cover, create subgraph using them
    from offsetbasedgraph.graphutils import GeneList
    genes = GeneList.from_file(args.genes_file_name).gene_list

    trans_regions = [g.transcription_region for g in genes]
    subgraph, trans = graph.create_subgraph_from_intervals(trans_regions, 20)

    # Translate genes using trans
    for g in genes:
        g.transcription_region = trans.translate(g.transcription_region)
        for exon in g.exons:
            exon = trans.translate(exon)

    levels = Graph.level_dict(subgraph.blocks)

    # Find start block by choosing a block having no edges in
    start = None
    for b in subgraph.blocks:
        if len(subgraph.reverse_adj_list[b]) == 0:
            start = b
            break

    print("== Subgraph ==")
    print(subgraph)

    assert start is not None

    from offsetbasedgraph import VisualizeHtml
    subgraph.start_block = start
    max_offset = sum([subgraph.blocks[b].length() for b in subgraph.blocks])
    v = VisualizeHtml(subgraph, 0, max_offset, 0, levels, "", 800, genes)
    print(v.get_wrapped_html())


def translate_genes_to_aligned_graph(args):
    trans = Translation.from_file(args.merged_graph_file_name)

    # Creat critical path multipath intervals of genes by translating
    # to complex graph.
    # Represent by txStart, txEnd as start,end and exons as critical intervals
    from offsetbasedgraph.graphutils import GeneList, translate_single_gene_to_aligned_graph
    genes = GeneList(get_gene_objects_as_intervals(args.genes, trans.graph1))
    gene_name = args.gene_name
    mpgenes = []
    mpgene_objects = []
    spgenes = []  # Also store single path for debugging
    n = 1
    for gene in genes.gene_list:
        n += 1
        if gene.name != gene_name:
            continue

        print("Found gene %s" % gene.name)
        print(gene)

        mpgene = translate_single_gene_to_aligned_graph(gene, trans)
        mpgene_objects.append(mpgene)
        #mpgenes.append(mpinterval)
        #mpgene_objects = MultiPathGene(gene.name, mpinterval)

        #spgenes.append(Gene(gene.name,
        #                    trans.translate(gene.transcription_region),
        #                    critical_intervals
        #                    )
        #               )


    import pickle

    mpgenes_list = GeneList(mpgene_objects)
    mpgenes_list.to_file(args.out_file_name)
    """
    with open("%s" % args.out_file_name, "wb") as f:
        pickle.dump(mpgenes, f)
    """

    for gene in spgenes:
        print("--------------")
        for exon in gene.exons:
            print(exon)

    gene_list = GeneList(spgenes)
    gene_list.to_file("%s_gene_list" % args.out_file_name)

    print(spgenes[0])
    print(spgenes[1])

    print("Genes written")


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
                #print("Match between")
                #print(g)
                #print(g2)
                equal += 1

            if g.faster_equal_critical_intervals(g2):
                #print("=== Exon match ===")
                #print(g)
                #print(g2)
                equal_exons += 1

    return equal, equal_exons

def analyse_multipath_genes2(args):
    import pickle

    from offsetbasedgraph.graphutils import create_gene_dicts, translate_single_gene_to_aligned_graph
    print("Reading genes")
    genes = GeneList(get_gene_objects_as_intervals(args.genes_file_name)).gene_list

    alt_loci_genes, gene_name_dict, main_genes = create_gene_dicts(genes)

    # alt loci genes are only genes on alt loci (nothing on main)
    # exon_dict contains only genes on main, index by offset of first exon

    # For every alt loci, create complex graph, translate genes and analyse them
    text_graph = create_initial_grch38_graph(args.chrom_sizes_file_name)
    graph, name_trans = grch38_graph_to_numeric(text_graph)

    orig_graph = graph.copy()
    name_trans_copy = name_trans.copy()

    equal_total = 0
    equal_exons_total = 0
    n_a = 1

    t1 = Translation.from_file("tmp_name_trans_single")
    t2 = Translation.from_file("name_trans_tmp")

    #assert t1 == name_trans

    n_alt_loci = len([bl for bl in text_graph.blocks if "alt" in bl])

    for b in text_graph.blocks:
    #for b in ["chr19_GL949747v2_alt", "chr19_KI270929v1_alt"]: #text_graph.blocks:
    #for b in ["chr13_KI270838v1_alt"]:
    #for b in ["chr19_KI270929v1_alt"]:
        if "alt" in b:
            print()
            #sys.stdout.write('\r  Analysing genes on alt locus ' + b + ' (number ' + str(n_a) + '/' + str(n_alt_loci) + ')' + ' ' * 20)
            print(b)
            sys.stdout.flush()

            n_a += 1
            genes_here = alt_loci_genes[b]
            trans, complex_graph = merge_alt_using_cigar(graph, name_trans, b)
            full_trans = name_trans + trans


            # Find candidates on main path to check against:
            genes_against = [g.copy() for g in main_genes[b]]
            genes_against_translated = []
            n = 0
            for mg in genes_against:
                #if n % 5 == 0 and n > 0:
                #    print("  Translated %d/%d genes" % (n, len(genes_against)))
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
            #print("\n  Candidates to check against: %d" % len(genes_against))

            """
            subgraph_intervals = []
            for g in genes_here_translated + genes_against_translated:
                subgraph_intervals.extend(g.critical_intervals)

            subgraph, trans, start_position = create_subgraph_around_alt_locus(complex_graph, full_trans, b, 30000)
            """

            equal, equal_exons = _analyse_multipath_genes_on_graph(genes_here_translated,
                                                                   genes_against_translated,
                                                                   complex_graph)
            equal_total += equal
            equal_exons_total += equal_exons

            #assert equal <= len(genes_here)
            #print("Total number of genes: %d" len(genes_here)), equal exons: %d, total %d genes" % (equal, equal_exons, len(genes_here)))

    print("SUM:")
    print("Equal: %d, equal exons: %d" % (equal_total, equal_exons_total))

    print("RESULTS:")
    print(" Number of genes originating from alt loci identical to a gene originaing from main chromosome: %d" % equal_total)
    print(" Number of genes with only identical exones (not start and end position): %d" % equal_exons_total)


def analyse_multipath_genes(args):

    import pickle
    with open("%s" % args.multipath_genes_file_name, "rb") as f:
        genes = pickle.loads(f.read())

    # Create a simple dict index to speed up duplicate search
    genes_index = {}
    for g in genes:
        if g.start_pos.offset in genes_index:
            genes_index[g.start_pos.offset].append(g)
        else:
            genes_index[g.start_pos.offset] = [g]

    print(genes)

    # Exon index
    print("Creating exon index")
    exon_index = {}
    for g in genes:
        first_exon = g.critical_intervals[0]
        index = "%s,%s" % (first_exon.region_paths[0], first_exon.start_position.offset)
        if index in exon_index:
            exon_index[index].append(g)
        else:
            exon_index[index] = [g]

    print("Created exon index")


    equal = 0
    equal_exons = 0
    n = 1
    for g in genes:
        if n % 1000 == 0:
            print("Checked %d genes" % n)
        n += 1
        #if g.start_pos.region_path_id != g.end_pos.region_path_id:
        #    print(g)

        #for g2 in genes_index[g.start_pos.offset]:

        first_exon = g.critical_intervals[0]
        index = "%s,%s" % (first_exon.region_paths[0], first_exon.start_position.offset)
        for g2 in exon_index[index]:#genes_index[g.start_pos.offset]:
            if g is g2:
                continue

            if g == g2:
                #print("=== Equal ==")
                #print(g)
                #print(g2)
                equal += 1

            if g.faster_equal_critical_intervals(g2):
                #print("== Equal critical intervals ==")
                #print(g)
                #print(g2)
                equal_exons += 1

    print("Equal: %d" % (equal / 2))
    print("Equal exons: %d" % (equal_exons / 2))


def html_alt_loci_select(args):
    # Prints all regions as an html select field
    html_out = """<select name='region'
               class='form-control' style='width: 320px;'>"""
    from offsetbasedgraph.graphutils import get_alt_loci_positions
    loci = get_alt_loci_positions()

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
