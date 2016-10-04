from offsetbasedgraph import LinearInterval, RegionPath, OffsetBasedGraph, GraphInterval

def make_lin_refs(mode, length, main_start, alt_start):
    lin_refs = {}
    if mode in ["identity", "poly", "delete"]:
        lin_refs["main"] = LinearInterval(
            "hg38",
            "main",
            main_start,
            main_start + length)
    if mode in ["identity", "poly", "insert"]:
        lin_refs["alt"] = LinearInterval(
            "hg38",
            "alt",
            alt_start,
            alt_start + length)

    return lin_refs


def test_graph(seq1, seq2, blocks, block_edges):
    seq1 = seq1.replace("-", "")
    seq2 = seq2.replace("-", "")
    cur_block = blocks[2]
    seqs = []
    name = "alt"
    for _ in range(10):
        cur_lin_ref = cur_block.linear_references[name]
        print (cur_lin_ref.start, cur_lin_ref.end)
        seqs.append(seq2[cur_lin_ref.start:cur_lin_ref.end])
        cur_block = blocks[[edge for edge in block_edges[cur_block.id]
                            if name in blocks[edge].linear_references][0]]
    seq = "".join(seqs)
    print "--------------------------------"

    print seq == seq2[:len(seq)]


def make_graph(seq1, seq2):
    # modes = ["poly", "insert", "delete", "identity"]
    block_id = 1
    cur_len = 0
    cur_mode = None
    alt_idx = 0
    main_idx = 0
    blocks = {}
    alt_block_ids = []
    main_block_ids = []
    main_block_id = 0
    alt_block_id = 0

    for c, c2 in zip(seq1, seq2):
        if c=="-" and c2=="-":
            raise Exception()
        if c == c2:
            new_mode = "identity"
        elif c == "-":
            new_mode = "insert"
        elif c2 == "-":
            new_mode = "delete"
        else:
            new_mode = "poly"

        if cur_mode is None:
            cur_mode = new_mode
        if new_mode == cur_mode:
            cur_len += 1
            continue

        lin_refs = make_lin_refs(cur_mode, cur_len, main_idx, alt_idx)
        if "alt" in lin_refs:
            alt_idx += lin_refs["alt"].length()
        if "main" in lin_refs:
            main_idx += lin_refs["main"].length()

        if cur_mode == "poly":
            cur_block_id = "main"+str(main_block_id)
            main_block_ids.append(cur_block_id)
            blocks[cur_block_id] = RegionPath(cur_block_id, {"main": lin_refs["main"]})
            main_block_id += 1
            cur_block_id = "alt"+str(alt_block_id)
            alt_block_ids.append(cur_block_id)
            blocks[cur_block_id] = RegionPath(cur_block_id, {"alt": lin_refs["alt"]})
            alt_block_id += 1
        else:
            if "main" in lin_refs:
                this_block_id = "main" + str(main_block_id)
                main_block_ids.append(this_block_id)
                main_block_id += 1
            else:
                this_block_id = "alt" + str(alt_block_id)
                alt_block_id += 1
            if "alt" in lin_refs:
                alt_block_ids.append(this_block_id)
            blocks[this_block_id] = RegionPath(this_block_id, lin_refs)

        cur_mode = new_mode
        cur_len = 1

    block_edges = {}
    for (out_id, in_id) in zip(alt_block_ids[:-1], alt_block_ids[1:]):
        block_edges[out_id] = [in_id]

    for (out_id, in_id) in zip(main_block_ids[:-1], main_block_ids[1:]):
        if out_id not in block_edges:
            block_edges[out_id] = []
        block_edges[out_id].append(in_id)

    graph = OffsetBasedGraph.create_graph_from_blocks(blocks, block_edges)
    return graph


def parse_blast(filename):
    f = open(filename, "r")
    text = f.read()
    lines = text.split("Query")[1:]
    query_lines = [line.split("\n")[0].strip().split() for line in lines]
    subject_lines = [line.split("\n")[2].strip().split() for line in lines]
    q_seqs = [line[1] if len(line)==3 else line[0] for line in query_lines]
    s_seqs = [line[2] if len(line)==4 else line[1] for line in subject_lines]
    q_seq = "".join(q_seqs)
    s_seq = "".join(s_seqs)
    return make_graph(q_seq, s_seq)

if __name__ == "__main__":
    graph = parse_blast("chr13alt_align.txt")
    graph_interval = GraphInterval.linear_segment_to_graph(graph, "main", 0, 1000)
    print graph_interval
    graph_interval = GraphInterval.linear_segment_to_graph(graph, "alt", 0, 1000)
    print graph_interval
