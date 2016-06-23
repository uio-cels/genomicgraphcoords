
from LinearInterval import LinearInterval
from RegionPath import RegionPath
from config import *
from DbWrapper import DbWrapper


class OffsetBasedGraph():
    """
    A simple graph structure, based on an offset-based coordinate system.
    Only region-paths are stored, not single bases.
    """
    def __str__(self):
        elements = [str(block.id) + ": " + str(block)
                    for block in self.blocks.values()]
        for key, val in self.block_edges.iteritems():
            elements.append("%s: %s" % (str(key), str(val)))
        return "\n".join(elements)

    def __init__(self, name):
        self.blocks = {}
        self.block_edges = {}
        self.cur_id = 0
        self.cur_id_counter = {}  # Number of times a base id of block
                                  # has been used
        self.db = DbWrapper()

    def get_blocks(self, lin_ref):
        return filter(
            lambda block: block.contains(lin_ref),
            self.blocks.values()
            )

    def get_intersecting_blocks(self, lin_ref):
        return filter(
            lambda block: block.intersects(lin_ref),
            self.blocks.values()
            )

    def get_block(self, lin_ref):
        filtered = self.get_blocks(lin_ref)
        if not filtered:
            return None
        return filtered[0]

    def merge_lin_refs(self,  org_lin_ref, new_lin_ref):
        org_block = self.get_block(org_lin_ref)
        if org_block is None:
            return

        pre_ref, post_ref = self.split_block(org_block, org_lin_ref)
        in_edges = self.get_previous_blocks(org_block.id)
        out_edges = []
        if org_block.id in self.block_edges:
            out_edges = self.block_edges[org_block.id]

        lin_refs = [pre_ref, post_ref, new_lin_ref, org_lin_ref]
        blocks = [self._add_block({"hg38": lr}) for lr in lin_refs]

        pre_block = blocks[0]
        post_block = blocks[1]
        new_block = blocks[2]
        org_block_1 = blocks[3]
        if pre_block.is_empty():
            for p_block in in_edges:
                self.add_block_edge(p_block, org_block_1.id)
                self.add_block_edge(p_block, new_block.id)
            del self.blocks[pre_block.id]
        else:
            for p_block in in_edges:
                self.add_block_edge(p_block, pre_block.id)
            self.add_block_edge(pre_block.id,  new_block.id)
            self.add_block_edge(pre_block.id,  org_block_1.id)

        if post_block.is_empty():
            self.block_edges[org_block_1.id] = out_edges[:]
            self.block_edges[new_block.id] = out_edges[:]
            del self.blocks[post_block.id]
        else:
            self.add_block_edge(new_block.id, post_block.id)
            self.add_block_edge(org_block_1.id, post_block.id)
            self.block_edges[post_block.id] = out_edges[:]

        # Remove org_block
        if org_block.id in self.block_edges:
            del self.block_edges[org_block.id]
        del self.blocks[org_block.id]

        for p_block in in_edges:
            self.block_edges[p_block].remove(org_block.id)

    def split_block(self, block, lin_ref):
        """
        Splits a block at the given linear reference. Returns two linear
        references, one before and one after lin_ref
        """
        current_lin_refs = filter(
            lambda lr: lr.genome_id == lin_ref.genome_id,
            block.linear_references.values())
        assert len(current_lin_refs) == 1
        current_lin_ref = current_lin_refs[0]
        splitted = []
        splitted.append(
            LinearInterval(lin_ref.genome_id,
                           lin_ref.chromosome,
                           current_lin_ref.start,
                           lin_ref.start,
                           lin_ref.strand)
        )
        splitted.append(
            LinearInterval(
                lin_ref.genome_id,
                lin_ref.chromosome,
                lin_ref.end,
                current_lin_ref.end,
                lin_ref.strand)
        )
        return splitted

    def merge_linear_segments(self, lin_seg_1, lin_seg_2):
        """
        Merges the graph where the two linear segments are
        """
        # First find the two blocks where each linear segment is
        block1 = self.get_block(lin_seg_1)
        block2 = self.get_block(lin_seg_2)
        if (block1 is None or block2 is None):

            if block2 is None:
                if DEBUG: print "_-_____"
                if DEBUG: print self.get_intersecting_blocks(lin_seg_2)
                if DEBUG: print "?????????????"
            if block1 is None:
                if DEBUG: print "+++++++++++++"
                if DEBUG: print self.get_intersecting_blocks(lin_seg_1)
                if DEBUG: print "?????????????"

            if DEBUG: print lin_seg_1, block1

            if DEBUG: print lin_seg_2, block2
            return
        # The new block that is a merge of two linear segments
        mid_block = self._add_block({lin_seg_1.genome_id: lin_seg_1,
                                     lin_seg_2.genome_id + "2": lin_seg_2})

        pre1, post1 = self.split_block(block1, lin_seg_1)
        pre2, post2 = self.split_block(block2, lin_seg_2)
        pre_block1 = RegionPath(block1.id, {pre1.genome_id: pre1})
        pre_block2 = RegionPath(block2.id, {pre2.genome_id: pre2})

        self.blocks[pre_block1.id] = pre_block1
        self.blocks[pre_block2.id] = pre_block2
        post_block1 = self._add_block({post1.genome_id: post1})
        post_block2 = self._add_block({post2.genome_id: post2})
        if block1.id in self.block_edges:
            self.block_edges[post_block1.id] = self.block_edges[block1.id][:]

        if block2.id in self.block_edges:
            self.block_edges[post_block2.id] = self.block_edges[block2.id][:]

        self.block_edges[pre_block1.id] = []
        self.block_edges[pre_block2.id] = []

        self.add_block_edge(pre_block1.id, mid_block.id)
        self.add_block_edge(pre_block2.id, mid_block.id)
        self.add_block_edge(mid_block.id, post_block1.id)
        self.add_block_edge(mid_block.id, post_block2.id)
        if DEBUG: print "-", lin_seg_1, lin_seg_2

    def add_block_edge(self, block_id, next_block_id):
        if block_id not in self.blocks:
            raise Exception("Edge id (%s) not in blocks: %s"
                            % (block_id, self.blocks))

        if block_id not in self.block_edges:
            self.block_edges[block_id] = [next_block_id]
            return

        self.block_edges[block_id].append(next_block_id)

    def get_subgraph(self, linear_interval, correct_alt):
        """Warning. Modifies current graph."""
        blocks = self.get_intersecting_blocks(linear_interval)
        min_C = 4000000000
        max_C = 0
        min_block = None
        max_block = None
        for block in blocks:
            lin_seg = filter(
                lambda li: li.chromosome == linear_interval.chromosome,
                block.linear_references.values())[0]
            if lin_seg.start < min_C:
                min_block = block
                min_C = lin_seg.start
            if lin_seg.end > max_C:
                max_block = block
                max_C = lin_seg.start

        subgraph = OffsetBasedGraph("pruned graph")

        cur_block = min_block
        subgraph.blocks[cur_block.id] = cur_block
        subgraph.start_block = min_block.id

        min_block.linear_references.values()[0].start = linear_interval.start
        max_block.linear_references.values()[0].end = linear_interval.end
        if min_block == max_block:
            return subgraph
        while True:
            found = False
            for new_block in self.block_edges[cur_block.id]:
                if any(["alt" in lr.chromosome and not lr.chromosome == correct_alt
                        for lr in self.blocks[new_block].linear_references.values()]):
                    continue
                subgraph.blocks[new_block] = self.blocks[new_block]
                if self.blocks[new_block] == max_block:
                    found = True
                    break
            if found:
                break
            if not self.block_edges[cur_block.id]:
                break
            cur_block = self.blocks[self.block_edges[cur_block.id][0]]

        for block in subgraph.blocks:
            subgraph.block_edges[block] = filter(
                lambda b: b in subgraph.blocks,
                self.block_edges[block])

        return subgraph

    def get_previous_blocks(self, block_id):
        previous_blocks = []
        for p_block, edges in self.block_edges.iteritems():
            if block_id in edges:
                previous_blocks.append(p_block)
        return previous_blocks

    def merge_from_chain_file(self, filename, genome1_id, genome2_id):
        segs = self._chain_file_to_linear_segments(
            filename, genome1_id, genome2_id)
        for seg1, seg2, in segs:
            self.merge_linear_segments(seg1, seg2)

    def add_species(self, genome_id, fn_chrom_sizes):
        """
        Adds blocks, one for each chromosome specified in the file
        """
        f = open(fn_chrom_sizes)
        for line in f.readlines():
            d = line.split()
            if "_" not in d[0]:
                meta = {}
                length = int(d[1])
                meta[genome_id] = LinearInterval(
                    genome_id, d[0], 0, length, "+")
                self._add_block(meta)

    def add_chromosome(self, genome_id, chromosome, chromosome_size):
        """
        Adds a whole chromosome as a linear block
        """
        linear_references = {}
        linear_references[genome_id] = LinearInterval(
            genome_id, chromosome, 0, chromosome_size, "+")
        self._add_block(linear_references)

    def _generate_id_from_linear_references(self, lrs):

        suggestion = ""
        sorted_genome_ids = sorted([lr for lr in lrs])
        if len(sorted_genome_ids) > 1:
            suggestion += "merge-"
        for genome_id in sorted_genome_ids:
            #suggestion += genome_id + "." + lrs[genome_id].chromosome + "-"
            # Skip non-alternative
            if len(sorted_genome_ids) == 1 or "alt" in lrs[genome_id].chromosome:
                suggestion += lrs[genome_id].chromosome + "-"

        base_id = suggestion

        if suggestion in self.cur_id_counter:
            suggestion += str(self.cur_id_counter[base_id] + 1)
            self.cur_id_counter[base_id] += 1
        else:
            suggestion += "0"
            self.cur_id_counter[base_id] = 0

        return suggestion

    def _add_block(self, linear_references):
        """
        Adds block for one genome
        linear_referes is a dict, having keys that are genome id and values
        that are LinearInterval objects
        """
        id = self._generate_id_from_linear_references(linear_references)
        block = RegionPath(id, linear_references)
        self.blocks[id] = block
        return block

    @classmethod
    def create_graph(cls, chrom_sizes, alt_loci_infos):
        """
        Creates and returns a block graph from two dicts
        """

        graph = cls("hg38")

        for chromosome, length in chrom_sizes.iteritems():
            graph.add_chromosome("hg38", chromosome, length)

        for info in alt_loci_infos:
            main_interval = LinearInterval(
                "hg38", info["chrom"], info["chromStart"], info["chromEnd"])

            alt_interval = LinearInterval(
                "hg38", info["name"], 0, info["length"])

            graph.merge_lin_refs(main_interval, alt_interval)

        return graph

    def include_alignments(self, alignments):
        for lr1, lr2 in alignments:
            if DEBUG: print "+", lr1, lr2
            self.merge_linear_segments(lr1, lr2)

    def pretty_alt_loci_name(self, id):
        """
        :return: Returns a prettier region name for an alternative loci
        """
        id = id.replace("hg38.", "")
        id = id.replace("merge-", "")
        if not "alt" in id:
            return id
        else:
            alt_id = id.split("-")[0]
            id = id.replace(alt_id, self.db.alt_loci_pretty_name(alt_id))
            #id += "test: ---%s---%s**" % (alt_id, self.db.alt_loci_pretty_name(alt_id))

        return id