"""
Methods for interacting with the different scripts. This file is called
via the web-server api, and prints html that is presented as results in the
tool.

run example:
$ python interface.py align_region2 chr7_KI270808v1_alt

"""
from __future__ import print_function
from __future__ import absolute_import

import sys
sys.path.append("/home/ivarandknut/python-projects/")
sys.path.append("/home/ivarandknut/python-projects/OffsetBasedGraph/")


import subprocess
import sys
import os
from config import *
from Visualize import Visualize, VisualizeHtml

from main import *

method = sys.argv[1]


if method == "test":
    import sys
    print(sys.version)
    print("Test is working .. ")


elif method == "get_all_regions":
    # Prints all regions as an html select field
    html_out = """<select name='region'
               class='form-control' style='width: 320px;'>"""
    db = DbWrapper()
    regions = db.get_alt_loci_infos()
    graph = OffsetBasedGraph("dummy-graph") # Only to access some methods

    for region in regions:
        if region["chromEnd"] - region["chromStart"] < 400000:
            html_out += "<option value='%s'>%s (%s:%d-%d)</option>" % \
                        (region["name"],
                         graph.pretty_alt_loci_name(region["name"]),
                         region["chrom"], \
                         region["chromStart"], region["chromEnd"])
    html_out += "</select>"
    print(html_out)


elif method == "align_region":
    """
    Aligns a given region, and visualizes using matplotlib.
    """
    db = DbWrapper()
    region = sys.argv[2]
    region_info = db.alt_loci_info(region)

    try:
        graph, orig_graph, gene_segments = create_align_graph(region, 0)
    except Exception as e:
        print(str(e)) # "<div class='alert alert-warning'>" + str(e) + "</div>"
        exit()

    area_start = region_info["chromStart"] - 50000
    area_end = region_info["chromEnd"] + 50000
    chrom = region_info["chrom"]
    graph = graph.get_subgraph(LinearInterval("hg38", chrom, area_start, area_end), region)

    description = "Visualization of graph created around <i>" + graph.pretty_alt_loci_name(region) + "</i>"
    #if len(gene_segments) > 3:
    #    gene_segments = gene_segments[0:3]

    v = Visualize(graph)
    v.show()


elif method == "align_region2" or method == "align_region_html":
    """
    Aligns a given region, and calls the Visualize class to print visualization
    of the results as html.
    """
    db = DbWrapper()
    region = sys.argv[2]
    region_info = db.alt_loci_info(region)

    #try:
    graph, orig_graph, gene_segments, gene_segments_orig = \
                                            create_align_graph(region, 0)
    #except Exception as e:
    #    print "<div class='alert alert-warning'>" + str(e) + "</div>"
    #    exit()

    area_start = region_info["chromStart"] - 50000
    area_end = region_info["chromEnd"] + 50000
    chrom = region_info["chrom"]
    graph = graph.get_subgraph(LinearInterval("hg38", chrom, area_start, area_end), region)



    description = "Visualization of graph created around <i>" + graph.pretty_alt_loci_name(region) + "</i>"
    #if len(gene_segments) > 3:
    #    gene_segments = gene_segments[0:3]

    """
    gene_segments2 = []
    for gene in gene_segments:
        name = gene.name
        g = db.get_gene(name)
        l = LinearInterval("hg38", g["chrom"], g["chromStart"], g["chromEnd"])
        s = Interval()
        segment = linear_segment_to_graph(graph, orig_graph, gene["chrom"],
                                              int(ex), int(ex_ends[ei]))
        s.create_from_block_list(graph.blocks[g["chrom"], int(g["chromStart"]), int(g["chromEnd"])])
        gene_segments2.append(s)

    genes = db.
    """



    v = VisualizeHtml(graph, area_start, area_end, 0, description, 800, gene_segments)
    print(str(v))
    print("<br><br>")

    total_width = v.width_used  #  Send the width used by the first visualizer
                                #  to the next

    subgraph_orig = orig_graph.get_subgraph(LinearInterval("hg38", chrom, area_start, area_end), region)
    v2 = VisualizeHtml(subgraph_orig, area_start, area_end, 1, "Original GRCh38 graph in the same area (no merges)", total_width + 60, gene_segments_orig)
    print(str(v2))

    # Print analysis details
    import globals
    print("<br><hr><br><h3>Details</h3>")
    print("""
    <p>The following is a brief description of how the topmost figure was created.
    The source code can be found <a href='https://github.com/uio-cels/gen-graph-coords' target='_blank'>here</a>.</p>
    <ol>
    """)
    print("""<li>Information about the alternative locus (region selected in the tool)
    was collected from the UCSC
    mysql database, using this command:<br>
    <i>mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
    -D hg38 -e "SELECT * FROM altLocations where name LIKE '%s' LIMIT 1;"</i>
    """ % region_info["name"])

    print("""<li>The sequence of the the alternative locus <i>%s</i> was fetched using
    togows.org api, by calling: <a href='%s' target='_blank'>%s</a></li>""" \
        % (region_info["name"], globals.togows_alt_url, globals.togows_alt_url))

    print("""<li>The sequence of the main path on GRCh38 corresponding to the alternative locus (%s) was fetched using
    the togows.org api, by calling: <a href='%s' target='_blank'>%s</a></li>""" \
        % (region_info["name"], globals.togows_main_url, globals.togows_main_url))

    from blast import BLAST_COMMAND
    blast_command = BLAST_COMMAND % ("main_path_sequence.fasta", "alt_locus_sequence.fasta", 95)
    print("""<li>The flanking regions in these two sequences were found as the sequence of base pairs at the beginning or end that were identical.</li>""")

    print(""" <li>The alternative locus was then merged with the main path where flanking regions were found.</li>
    <li>Genes were collected from the UCSC hg38
    database, table: <i>knownGene</i>.
    Only genes present in multiple region paths in the part of the graph visualized were selected.</li>
    """)
    print("</ol>")