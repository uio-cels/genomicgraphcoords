"""
A collections of small scripts used to do various tasks by calling the classes.
"""
from blast import blast_align, get_filtered_alignments
from TogowsWrapper import save_sequence_to_fasta
from OffsetBasedGraph import OffsetBasedGraph
from DbWrapper import DbWrapper
from LinearInterval import LinearInterval
from Visualize import Visualize
from config import *
from os.path import isfile
from CoordinateConverter import create_block_index, linear_segment_to_graph
from Interval import Interval, IntervalCollection


def get_intereseting_points_from_alignments(alignments):
    if not alignments:
        return {}
    interresting_points = {alignments[0][0].chromosome: [],
                           alignments[0][1].chromosome: []}
    for alignment in alignments:
        interresting_points[alignment[1].chromosome].append(
            (alignment[1].start, alignment[1].end))
        interresting_points[alignment[0].chromosome].append(
            (alignment[0].start, alignment[0].end))
    return interresting_points


def get_alignments(region_name, alt_info):
    align_file = DATA_PATH + "alignment%s.tmp" % region_name
    import globals
    globals.blast_result = "alignment%s.tmp" % region_name
    alt_fasta = save_sequence_to_fasta(
        region_name, 0, alt_info["length"])
    consensus_fasta = save_sequence_to_fasta(
        alt_info["chrom"], alt_info["chromStart"], alt_info["chromEnd"])

    if isfile(align_file):
        return align_file  # File is cached already

    align_file = blast_align(consensus_fasta, alt_fasta, 95, 300, align_file)

    return align_file


def get_flanking_alignments(region_name, alt_info):
    align_file = DATA_PATH + "alignment%s.tmp" % region_name
    import globals
    globals.blast_result = "alignment%s.tmp" % region_name
    alt_fasta = save_sequence_to_fasta(
        region_name, 0, alt_info["length"])
    consensus_fasta = save_sequence_to_fasta(
        alt_info["chrom"], alt_info["chromStart"], alt_info["chromEnd"])

    alt_seq = "".join(map(lambda x: x.strip(), open(alt_fasta, "r").readlines()[1:]))
    consensus_seq = "".join(map(lambda x: x.strip(), open(consensus_fasta, "r").readlines()[1:]))
    start_flank_length = 0
    for i in xrange(min(len(consensus_seq), len(alt_seq))):
        if alt_seq[i] != consensus_seq[i]:
            start_flank_length = i
            break
    else:
        raise Exception("Main and Alt is completely equal")

    linRef1 = LinearInterval("hg38", alt_info["chrom"], alt_info["chromStart"],
                             alt_info["chromStart"] + start_flank_length)
    linRef2 = LinearInterval("hg38", region_name, 0, start_flank_length)
    start_pair = (linRef1, linRef2)
    stop_flank_length = 0
    print consensus_seq[-100:]
    print alt_seq[-100:]
    for i in xrange(min(len(consensus_seq), len(alt_seq))):
        if alt_seq[-i-1] != consensus_seq[-i-1]:
            stop_flank_length = i
            break
    else:
        raise Exception("Main and Alt is completely equal")
    print stop_flank_length
    stopRef1 = LinearInterval(
        "hg38", alt_info["chrom"],
        alt_info["chromEnd"]-stop_flank_length, alt_info["chromEnd"])
    stopRef2 = LinearInterval(
        "hg38", region_name,
        alt_info["length"]-stop_flank_length, alt_info["length"])
    stop_pair = (stopRef1, stopRef2)
    return [start_pair, stop_pair]


segments = []


def create_align_graph(region_name, min_length):
    dbw = DbWrapper()
    chromosome_ids = ["chr%s" % n for n in range(1, 23)+["X", "Y"]]
    chromosome_sizes = [dbw.chrom_length(chr_id) for chr_id in chromosome_ids]
    chromosome_info = dict(zip(chromosome_ids, chromosome_sizes))

    alt_infos = dbw.get_alt_loci_infos()

    graph = OffsetBasedGraph.create_graph(
        chromosome_info, alt_infos)

    import copy
    orig_graph = copy.deepcopy(graph)

    alt_info = dbw.alt_loci_info(region_name)




    #align_file = get_alignments(region_name, alt_info)
    #alignments = get_filtered_alignments(align_file, alt_info["chromStart"])
    alignments = get_flanking_alignments(region_name, alt_info)


    #print "Alignments :"
    #print alignments
    if len(alignments) == 0:
        raise Exception("There were no alignments between the alternative locus and consensus path in the the given region. Try another region.")

    if DEBUG: print "ALIGNMENTS"
    for a in alignments:
        if DEBUG: print a
    graph.include_alignments(alignments)

    interesting_points = get_intereseting_points_from_alignments(alignments)
    genes = []
    values = interesting_points.keys()
    for i in xrange(len(interesting_points[values[0]])):
        pgs0 = [dbw.genes_crossing_position(
            v, interesting_points[v][i][0]) for v in values]

        if all(pgs0):
            genes = [gs[0] for gs in pgs0]
            break

        pgs1 = [dbw.genes_crossing_position(
            v, interesting_points[v][i][1]) for v in values]
        if all(pgs1):
            genes = [gs[0] for gs in pgs1]
            break
    gene_intervals = {}
    exons = []
    block_graph = create_block_index(graph)
    block_graph_orig = create_block_index(orig_graph)
    for gene in genes:
        #print "Gene %s, %s" % (gene["gname"], gene["name"])
        l = LinearInterval(
            "hg38", gene["chrom"], gene["txStart"],
            gene["txEnd"], gene["strand"])
        gene_intervals[gene["name"]] = l
        l.gene_name = gene["gname"]

        #print gene["exonStarts"]
        #print gene ["exonEnds"]

        # Create one interval for each exon
        ei = 0
        ex_ends = gene["exonEnds"].split(",")
        for ex in gene["exonStarts"].split(","):
            if ex == "":
                continue

            #print "<p>Ex: %d, %d</p>" % (int(ex), int(ex_ends[ei]))
            l = LinearInterval("hg39", gene["chrom"], int(ex), int(ex_ends[ei]), "+")
            segment = linear_segment_to_graph(graph, block_graph, gene["chrom"],
                                              int(ex), int(ex_ends[ei]))
            segment.name = "Exon for " + gene["name"]
            segment.gene_name = gene["name"]
            segment.is_exon = True
            #print "<p>Exon segment: <p>"
            #print segment

            exons.append(segment)
            ei += 1



    for gi in gene_intervals.values():
        if DEBUG: print gi
        if DEBUG: print "###", graph.get_intersecting_blocks(gi)

    gene_segments = []
    gene_segments_orig_graph = []

    for name, interval in gene_intervals.iteritems():
        segment = linear_segment_to_graph(
                graph, block_graph,
                interval.chromosome, interval.start, interval.end)
        segment.name = name
        segment.gene_name = interval.gene_name
        if DEBUG: print segment
        segments.append(segment)
        gene_segments.append(segment)

        # Segment on original graph
        segment_orig = linear_segment_to_graph(
                orig_graph, block_graph_orig,
                interval.chromosome, interval.start, interval.end)
        segment_orig.name = name
        segment_orig.gene_name = interval.gene_name
        gene_segments_orig_graph.append(segment_orig)


    gene_segments = list(reversed(gene_segments[0:3])) # Limit to max three genes
    gene_segments_orig_graph = list(reversed(gene_segments_orig_graph[0:3])) # Limit to max three genes

    gene_segments.extend(exons)

    return graph, orig_graph, gene_segments, gene_segments_orig_graph


def is_good_loci(region_name, dbw):
    alt_info = dbw.alt_loci_info(region_name)
    align_file = get_alignments(region_name, alt_info)

    alignments = get_filtered_alignments(align_file, alt_info["chromStart"])
    interesting_points = get_intereseting_points_from_alignments(alignments)

    values = interesting_points.keys()
    for i in xrange(len(interesting_points[values[0]])):
        if all([dbw.genes_crossing_position(v, interesting_points[v][i][0])
                for v in values]):
            return True
        if all([dbw.genes_crossing_position(v, interesting_points[v][i][1])
                for v in values]):
            return True

    return False


def find_good_loci():
    dbw = DbWrapper()
    alt_infos = dbw.get_alt_loci_infos()
    names = [alt_info["name"] for alt_info in alt_infos]
    for name in names:
        if is_good_loci(name, dbw):
            if DEBUG: print name
    return

if __name__ == "__main__":

    graph = create_align_graph("chr11_KI270927v1_alt", 0)
    # chr11_KI270827v1_alt", 0)

    graph = graph.get_subgraph(
        LinearInterval("hg38", "chr11", 900000, 1300000),
        "chr11_KI270927v1_alt")


    if DEBUG: print "BLOCKS: "
    for b in graph.blocks:
        if DEBUG: print "BLock " + str(b) + " : " + str(graph.blocks[b])



    if DEBUG: print "EDGES: "
    for k, v in graph.block_edges.iteritems():
        if DEBUG: print k,"---", v

    for edge in graph.block_edges:
        assert(len(graph.block_edges[edge]) <= 2)

    v = Visualize(graph)
    v.visualize_intervals(segments)
    v.show()
#
#    interval = Interval(graph)
#    block_list = ["hg38.chr11-31", "hg38.chr11_KI270827v1_alt-0"]
#    interval.create_from_block_list(block_list, 10000, 50000)
#    v.visualize_intervals([interval])
#
#    block_list = ["hg38.chr11_KI270827v1_alt-0", "hg38.chr11-34"]
#    interval.create_from_block_list(block_list, 40000, 20000)
#    v.visualize_intervals([interval])
#

