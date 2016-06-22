"""
Methods for interacting with the different scripts.
"""
import subprocess
import sys
import os
import config
from Visualize import Visualize, VisualizeHtml

from main import *

method = sys.argv[1]

if method == "test":
    print "Test is working .. "
    print os.path.dirname(os.path.realpath(__file__))
    print DATA_PATH

    dir = "/var/www/web/data/tmp/"
    blast = 'blastn -outfmt "6 qseqid sseqid qstart qend sstart send length score bitscore evalue" -query %schr11:937258-1150312.fasta -subject %schr11_KI270927v1_alt:1-218613.fasta -perc_identity 95.00 ' % (dir, dir)
    p = subprocess.Popen(blast, shell=True,
                         cwd=os.path.dirname(os.path.realpath(__file__)),
                         stdout=subprocess.PIPE)
    p.wait()
    result = p.communicate()[0]
    print blast

    print result

elif method == "align_region":

    #graph = create_align_graph("chr11_KI270927v1_alt", 0)
    graph = create_align_graph("chr2_KI270894v1_alt", 0)
    area_start = 90000150
    area_end = 90450000
    #graph = graph.get_subgraph(LinearInterval("hg38", "chr11", 900000, 1200000), "chr11_KI270927v1_alt")
    #graph = graph.get_subgraph(LinearInterval("hg38", "chr2", area_start, area_end), "chr2_KI270894v1_alt")

    """
    print "BLOCKS: "
    for b in graph.blocks:
        print "BLock " + str(b) + " : " + str(graph.blocks[b])



    print "EDGES: "
    for k, v in graph.block_edges.iteritems():
        print k,"---", v

    for edge in graph.block_edges:
        assert(len(graph.block_edges[edge]) <= 2)
    """

    """
    v = Visualize(graph)
    v.visualize_intervals(segments)
    v.show()
    """
    v = VisualizeHtml(graph, area_start, area_end)
    print str(v)

elif method == "get_all_regions":
    html_out = "<select name='region' class='form-control' style='width: 320px;'>"
    db = DbWrapper()
    regions = db.get_alt_loci_infos()
    graph = OffsetBasedGraph("dummy-graph") # Only to access some methods

    for region in regions:
        if region["chromEnd"] - region["chromStart"] < 400000:
            html_out += "<option value='%s'>%s (%s:%d-%d)</option>" % \
                        (region["name"], graph.pretty_alt_loci_name(region["name"]), region["chrom"], \
                         region["chromStart"], region["chromEnd"])
    html_out += "</select>"

    print html_out


elif method == "align_region2":
    """
    Align region given by input parameter
    """
    db = DbWrapper()
    region = sys.argv[2]

    region_info = db.alt_loci_info(region)


    #graph = create_align_graph("chr11_KI270927v1_alt", 0)
    try:
        graph, orig_graph, gene_segments = create_align_graph(region, 0)
    except Exception as e:
        print "<div class='alert alert-warning'>" + str(e) + "</div>"
        exit()

    area_start = region_info["chromStart"] - 50000
    area_end = region_info["chromEnd"] + 50000
    chrom = region_info["chrom"]
    #graph = graph.get_subgraph(LinearInterval("hg38", "chr11", 900000, 1200000), "chr11_KI270927v1_alt")
    graph = graph.get_subgraph(LinearInterval("hg38", chrom, area_start, area_end), region)

    """
    print "BLOCKS: "
    for b in graph.blocks:
        print "BLock " + str(b) + " : " + str(graph.blocks[b])

    print "EDGES: "
    for k, v in graph.block_edges.iteritems():
        print k,"---", v

    for edge in graph.block_edges:
        assert(len(graph.block_edges[edge]) <= 2)
    """

    """
    v = Visualize(graph)
    v.visualize_intervals(segments)
    v.show()
    """
    description = "Visualization of graph created around <i>" + graph.pretty_alt_loci_name(region) + "</i>"
    if len(gene_segments) > 3:
        gene_segments = gene_segments[0:3]


    gene_segments2 = []
    for gene in gene_segments:
        name = gene.name
        g = db.get_gene(name)
        l = LinearInterval("hg38", g["chrom"], g["chromStart"], g["chromEnd"])
        s = Interval()
        s.create_from_block_list(graph.blocks[g["chrom"], int(g["chromStart"]), int(g["chromEnd"])])
        gene_segments2.append(s)

    genes = db.

    v = VisualizeHtml(graph, area_start, area_end, 0, description, 800, gene_segments)

    #print "<p>There are %d example genes in this area</p>" % len(gene_segments)
    #v.visualize_intervals(gene_segments[0:1])

    print str(v)
    print "<br><br>"

    total_width = v.width_used

    subgraph_orig = orig_graph.get_subgraph(LinearInterval("hg38", chrom, area_start, area_end), region)
    v2 = VisualizeHtml(subgraph_orig, area_start, area_end, 1, "Original GRCh38 graph", total_width + 60)
    print str(v2)

    import globals
    print "<br><hr><br><h3>Details</h3>"
    print """
    <p>The following is a brief description of how the topmost figure was created.
    The source code can be found <a href=''>here</a>.</p>
    <ol>
    """
    print """<li>Information about the alternative locus (region selected in the tool)
    was collected from the UCSC
    mysql database, using this command:<br>
    <i>mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
    -D hg38 -e "SELECT * FROM altLocations where name LIKE '%s' LIMIT 1;"</i>
    """ % region_info["name"]

    print """<li>The sequence of the this alternative locus <i>%s</i> was fetched using
    togows.org api, by calling: <a href='%s' target='_blank'>%s</a></li>""" \
        % (region_info["name"], globals.togows_alt_url, globals.togows_alt_url)

    print """<li>The sequence of the main path on GRCh38 corresponding to alternative locus <i>%s</i> was fetched using
    togows.org api, by calling: <a href='%s' target='_blank'>%s</a></li>""" \
        % (region_info["name"], globals.togows_main_url, globals.togows_main_url)

    from blast import BLAST_COMMAND
    blast_command = BLAST_COMMAND % ("main_path_sequence.fasta", "alt_locus_sequence.fasta", 95)
    print """<li>These two sequences were aligned using blastn.
    The following command was run:<br>%s.
    <br>The result <a href='http://46.101.93.163/data/tmp/%s' target='_blank'>
    can be seen here</a>.</li>""" % (blast_command, globals.blast_result)
    print """<li>For this example, only alignments in which both sequences
    had 50 or more bp were chosen. Alignments were chosen from the top
    of the file, excluding alignments that were either overlapping or
    <i>crossing</i> previously selected alignments. With <i>crossing</i> we mean
    an aligment between two sequences, where one of the two sequences
    is before and the other after a sequence aligned previously.</li>
    <li>The alternative locus was then merged with the main path based on
     these alignments.</li>
    <li>Genes (if any were visualized) were collected from the UCSC hg38
    database, table: <i>knownGene</i>. Only genes that were crossing
    'bundaries' in the graph where selected. With crossing a boundary,
    we mean genes that had start position before a visualized edge
    and end position after a visualized edge in the graph.</li>
    """




    print "</ol>"