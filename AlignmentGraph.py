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
    out_node = int(parts[1])
    in_node = int(parts[3])
    return (out_node, in_node)


def parse_gfa(filename):
    graph = OffsetBasedGraph("alignment-graph")
    f = open(filename, "r")
    cur_id = 0
    cur_len = None
    cur_intervals = []
    for line in f.readlines():
        parts = line.split()
        if parts[0] == "H":
            continue
        _id = int(parts[1])
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
    graph.blocks[cur_id] = RegionPath(
        cur_id, {li.chromosome: li for li in cur_intervals})

    return graph


def update_offsets(graph):
    start_blocks_ids = [block_id for block_id in graph.blocks if
                        not any([block_id in in_edges for
                                 in_edges in graph.block_edges.values()])]
    for start_block_id in start_blocks_ids:
        print(start_block_id)
        chromosomes = graph.blocks[start_block_id].linear_references.keys()
        for chromosome in chromosomes:
            print(chromosome)
            cur_block = graph.blocks[start_block_id]
            offset = 0
            while True:
                offset += cur_block.get_length()
                if cur_block.id not in graph.block_edges:
                    break
                next_blocks = [
                    e for e in graph.block_edges[cur_block.id] if
                    chromosome in graph.blocks[e].linear_references]
                if not next_blocks:
                    break
                cur_block = graph.blocks[next_blocks[0]]
                lin_ref = cur_block.linear_references[chromosome]
                lin_ref.start += offset
                lin_ref.end += offset


if __name__ == "__main__":
    graph = parse_gfa("tmp.gfa")
    for v in graph.blocks.values():
        print(v)
    print(graph.block_edges)
    update_offsets(graph)
    for v in graph.blocks.values():
        print(v)
