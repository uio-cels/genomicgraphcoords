
# Live demo

For a live demo of this tool, please follow this link:  http://46.101.93.163/gen-graph-coords/


# About
A tool to merge segments of alternative loci with the main paht of GRCh38 and find genes that covers both merged and unmerged paths.
class OffsetBasedGraph is a simple implementation of a sequence graph with offset based coordinates
class LinearInterval represents intervals on linear reference genomes
class Interval represents intervals on graph based reference genomes, as described in the article
class RegionPath represents region paths in the graph, as described in the article

 To run from commandline:
> python interface.py chrN_GENBANKSEQACCS_alt
Example:
> python interface.py chr7_KI270808v1_alt


...