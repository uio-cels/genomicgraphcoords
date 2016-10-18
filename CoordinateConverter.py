"""
Various methods for converting linear coordinates to a graph
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from Interval import Interval, IntervalCollection
from config import *

def create_block_index(graph):
    """ Returns a dictionary where indices are linear block ids (chromosome ids
    and alt loci ids), and values are list of blocks in the graph that come
    from these blocks.
    """

    index = {}
    for block_id, block in graph.blocks.items():
        for linear_reference in list(block.linear_references.values()):
            chr_id = linear_reference.chromosome
            if chr_id in index:
                index[chr_id].append(block)
            else:
                index[chr_id] = [block]

    return index


def linear_segment_to_graph(graph, graph_block_index, chr_id, start, end):
    """
    Takes a linear segment on hg38 and returns a ChainSegment object
    """
    if DEBUG: print("Linear segment to graph %s, %d, %d" % (chr_id, start, end))

    block_list = []
    start_pos = 0
    end_pos = 0

    start_block, start_pos = linear_coordinate_to_graph(graph, graph_block_index, chr_id, start)
    end_block, end_pos = linear_coordinate_to_graph(graph, graph_block_index, chr_id, end)
    block_list.append(start_block)
    current_block = start_block
    while True:
        # Follow edge to next block. If there are more than one edge, follow
        # edge to the block that is not an alternative loci
        if current_block == end_block:
            if DEBUG: print("Current block == end block")
            break
        prev_block = current_block
        edges = graph.block_edges[current_block]

        # Find next current block. Either alt of consensus
        for edge in edges:

            if DEBUG: print("edge" + str(edge))
            if "alt" in chr_id and "alt" in edge:
                current_block = edge
                if DEBUG: print("    GOing to alt")
                break
            elif not "alt" in chr_id and not "alt" in edge:
                current_block = edge
                break
            else:
                current_block = edge

        if current_block == prev_block:
            raise Exception("Error while traversing block. Did not find " +\
                            "next block for block %s, %s" % (current_block,edges))
        block_list.append(current_block)


    seg = Interval(graph)
    seg.create_from_block_list(block_list, start_pos, end_pos)
    return seg


def linear_coordinate_to_graph(graph, graph_block_index, chr_id, coordinate):
    """ Maps the linear coordinate to a coordinate in the block graph.
    graph_block_index should be sent as a parameter, and can be obtained
    by calling create_block_index()
    """
    if DEBUG: print("coordinate to graph. " + str(chr_id) + " " + str(coordinate))
    # Our coordinate can be in any of these blocks
    print("block index")
    print(graph_block_index.keys())
    potential_blocks = graph_block_index[chr_id]

    for potential_block in potential_blocks:
        # Get start and end position of this block in linear genome,
        # check whether this is the correct block
        # assume always only one linear refernce, i.e. only one species
        linear_references = list(potential_block.linear_references.values())
        for lr in linear_references:

            start = lr.start
            end = lr.end
            #if DEBUG: print "     %s, %d, %d " % (lr.chromosome, start, end)

            if start <= coordinate and end >= coordinate:
                return (potential_block.id, coordinate - start)

    raise Exception("No block found for chr_id %s, coordinate %d" % (chr_id, coordinate))

def graph_to_linear():
    pass
