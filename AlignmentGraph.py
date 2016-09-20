from LinearInterval import LinearInterval
from RegionPath import RegionPath
from OffsetBasedGraph import OffsetBasedGraph


def parse_sequence(parts):
    seq = parts[2]
    return len(seq)


def parse_path(parts, l):
    chromosome = parts[2]
    strand = parts[4]
    return LinearInterval("hg38", chromosome, 0, l, strand)


def parse_edge(parts):
    out_node = parts[1]
    in_node = parts[3]
    return (out_node, in_node)


def parse_gfa(filename):
    graph = OffsetBasedGraph("alignment-graph")
    f = open(filename, "r")
    cur_id = 0
    cur_len = None
    cur_intervals = []
    for line in f.readlines():
        parts = line.split()
        _id = parts[1]
        if cur_id != _id:
            if cur_len is not None:
                graph.blocks[cur_id] = RegionPath(
                    cur_id, {li.chromosome: li for li in cur_intervals})
            cur_len = None
            cur_intervals = []
            cur_id = _id
        if parts[0] == "S":
            cur_len = parse_sequence(parts)
        elif parts[0] == "P":
            cur_intervals.append(parse_path(parts, cur_len))
        elif parts[0] == "L":
            out_id, in_id = parse_edge(parts)
            if out_id not in graph.block_edges:
                graph.block_edges[out_id] = []
            graph.block_edges[out_id].append(in_id)

    return graph

if __name__ == "__main__":
    graph = parse_gfa("tmp.gfa")
    print graph.blocks.keys()
    print graph.block_edges
