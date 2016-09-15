from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import numpy as np

class Interval(object):
    """
    Class for genomic interval on a graph-based reference genome

    """

    def __init__(self, rg):
        """
        Inits an empty interval
        :param rg: The reference graph this interval is on. Typically
        an OffsetBasedGraph object.
        """
        self.reference_graph = rg
        self.block_list = None
        self.start_block = None
        self.end_block = None
        self.start_pos = None
        self.end_pos = None
        self.meta = ""
        self.name = ""
        self.gene_name = ""
        self.is_exon = False



    def create_from_block_list(self, block_list, start_pos, end_pos):
        """
        Builds the interval from a list of region path identifiers, and a start
        and end position.
        :param block_list: List of region path identifiers
        :param start_pos:
        :param end_pos:
        """
        self.block_list = block_list
        self.start_block = block_list[0]
        self.end_block = block_list[-1]
        self.start_pos = start_pos
        self.end_pos = end_pos

    def overlap(self, other):
        """
        Computes the number of base pairs overlap with another interal
        """
        overlap = 0
        for block in [self.reference_graph.blocks[block_id] for block_id in self.block_list]:
            if block not in [self.reference_graph.blocks[block_id] for block_id in other.block_list]:
                continue
            start = 0
            end = block.length()
            if block == self.reference_graph.blocks[self.start_block]:
                start = max(self.start_pos, start)
            if block == self.reference_graph.blocks[other.start_block]:
                start = max(other.start_pos, start)
            if block == self.reference_graph.blocks[self.end_block]:
                end = min(self.end_pos, end)
            if block == self.reference_graph.blocks[other.end_block]:
                end = min(other.end_pos, end)
            overlap += end-start
        return overlap

    def contains_block(self, block):
        return block in self.block_list

    def __eq__(self, other):
        """
        :param other: An other interval
        :return: Returns true if this interval is the same as the other,
        meaning it contains the exact same base pairs.
        """
        if self.start_pos != other.start_pos:
            return False
        if self.end_pos != other.end_pos:
            return False
        if self.start_block != other.start_block \
                or self.end_block != other.end_block:
            return False
        if len(self.block_list) != len(other.block_list):
            return False
        for i in range(0, len(self.block_list)):
            if self.block_list[i] != other.block_list[i]:
                return False
        return True

    def create_from_file_line(self, line):
        """
        Builds the inteval from a line in a a interval collection file

        Assumed convention for a interval line:
        first_block_offset first block second_block,...,last_block_offset
            last_block_end_pos  metadata1 metadata2

        Blocks are separated by spaces, and metadata by tabs

        :param line: A line (string) from a interval collection file
        :return: None
        """

        self.block_list = []
        line = line.split("\t")
        ids = line[0].split(" ")
        self.start_block = ids[1]
        self.start_pos = int(ids[0])
        for inner_block in ids[1:-1]:
            self.block_list.append(inner_block)
        self.end_block = self.block_list[-1]
        self.end_pos = int(ids[-1])

    def to_file_line(self):
        line = str(self.start_pos)  + " "
        for block in self.block_list:
            line += block + " "
        line += str(self.end_pos)
        line += "\t+\t100\tmetadata"
        return line

    def get_ids(self):
        ids = []
        for block in self.block_list:
            start = 0
            end = self.reference_graph.blocks[block].length()
            if block == self.start_block:
                start = self.start_pos
            if block == self.end_block:
                end = self.end_pos
            print("Start" , start)
            print(end)
            ids += [block + "-" + str(i) for i \
                    in range(start, end)]

        return ids

    def __str__(self):
        s = str(self.start_pos)
        for block in self.block_list:
            s += ", " + self.reference_graph.pretty_alt_loci_name(block)
        s += ", " + str(self.end_pos)

        return s
        #return str(self.start_pos) + " " + ", ".join(self.block_list) + ",  " \
        #            + str(self.end_pos)

    def __repr__(self):
        return self.__str__()


class IntervalCollection(object):
    def __init__(self, reference_graph):
        self.reference_graph = reference_graph
        self.segments = []

    def create_from_file(self, filename):
        f = open(filename)
        for line in f.readlines():
            if line[0] != "#":  # Allowing comment lines starting with #
                seg = Interval(self.reference_graph)
                seg.create_from_file_line(line)
                self.segments.append(seg)

    def write_to_file(self, filename):
        f = open(filename, "w")
        lines = []
        for segment in self.segments:
            lines.append(segment.to_file_line() + "\n")
        f.writelines(lines)

    def __len__(self):
        return len(self.segments)

    def __getitem__(self, key):
        return self.segments[key]

    def create_from_segment_list(self, segment_list):
        self.segments = segment_list

    def overlap(self, other):
        """
        Computes the overlap with another segment collection
        :param other:
        :return: The number of bases that overlap between the two collections
        """
        overlap = 0
        for seg1 in self.segments:
            for seg2 in other.segments:
                overlap += seg1.overlap(seg2)
        return overlap

    def size(self):
        n_bases = 0
        for seg in self.segments:
            n_bases += seg.size()
        return n_bases

    def get_ids(self):
        ids = []
        for seg in self.segments:
            ids += seg.get_ids()
        return ids

