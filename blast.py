import subprocess
import os
from config import *
from LinearInterval import LinearInterval

BLAST_COMMAND = 'blastn -outfmt "6 qseqid sseqid qstart qend sstart send length score bitscore evalue" -query %s -subject %s -perc_identity %.2f -num_threads 100'


class Alignment:
    """
    Simple class that represents one alignment
    """
    def __init__(self, chr1, chr2, start1, end1, start2, end2, score):
        self.chr1 = chr1
        self.chr2 = chr2
        self.start1 = start1-1
        self.start2 = start2-1
        self.end1 = end1-1
        self.end2 = end2-1
        self.score = score

    def __repr__(self):
        return "%s:%d-%d --- %s:%d:%d\n" % (self.chr1, self.start1, self.end1,
                                            self.chr2, self.start2, self.end2)


def blast_align(fasta1, fasta2, min_ratio_match, min_length,
                out_name="alignment.tmp"):

    if DEBUG: print "___ RUNNING BLAST_ALIGN _____ "

    blast_command = BLAST_COMMAND % (fasta1, fasta2, min_ratio_match)
    if DEBUG: print "\n " + blast_command + " \n"
    p = subprocess.Popen(blast_command, shell=True,
                         cwd=os.path.dirname(os.path.realpath(__file__)),
                         stdout=subprocess.PIPE)


    #if DEBUG: print "Path 1: " + str(os.path.abspath(os.curdir))

    p.wait()
    result = p.communicate()[0]

    f = open(out_name, "w")
    f.write(result)
    f.close()

    return out_name


def _is_crossing(a1, a2):
    """ Returns true if the alignment a1 is crossing a2
    """
    if a1.chr1 != a1.chr1:
        return False
    if a1.start1 < a2.start1 and a1.start2 > a2.start2:
        return True
    if a2.start1 < a1.start1 and a2.start2 > a1.start2:
        return True


def _is_overlapping(a1, a2):
    """
    Return true if alignment a1 is overlapping a2
    """
    if a1.chr1 != a2.chr1:
        return False

    if (a1.start1 < a2.start1 and a1.end1 > a2.start1) or (a1.end1 > a2.end1 and a1.start1 < a2.end1):
        return True

    if (a1.start2 < a2.start2 and a1.end2 > a2.start2) or (a1.end2 > a2.end2 and a1.start2 < a2.end2):
        return True


def _alignment_is_valid(alignments, a):
    """
    Returns true if alignemtn a is valid according to alignment
    """
    for a2 in alignments:
        if _is_crossing(a2, a) or _is_overlapping(a2, a):
            return False
    return True


def filter_alignments(alignments, min_length=50, non_overlapping=True,
                      no_crossing=True):
    """
    :param min_length:
    :param non_overlapping:
    :param no_crossing:
    :return:
    """

    # Alignments are already sorted from best to worst.
    # Keep all, except those who are reverse, too short, are overlapping
    filtered = []
    for a in alignments:
        if a.end1 - a.start1 < min_length or a.end2 - a.start2 < min_length:
            break
        if not _alignment_is_valid(filtered, a):
            break

        # Temporarily: Do not allow allignements at beginning or end
        #if a.start1 == 0 or a.start2 == 0:
        #    break
        filtered.append(a)


    return filtered


def get_alignments(filename, offset=0):
    alignments = []

    f = open(filename, "r")
    print "Get alignments"
    for line in f.readlines():
        print "line"
        parts = line.split()
        chrId1 = parts[0].split(":")[1]
        chrId2 = parts[1].split(":")[1]
        strt1 = int(parts[2])+offset
        end1 = int(parts[3])+offset
        strt2 = int(parts[4])
        end2 = int(parts[5])

        a = Alignment(chrId1, chrId2, strt1, end1, strt2, end2, int(parts[6]))
        alignments.append(a)
        f.close()

    return alignments


def create_linear_intervals_from_alignments(alignments):
    lin_ref_pairs = []
    for a in alignments:
        linRef1 = LinearInterval("hg38", a.chr1, a.start1, a.end1)
        linRef2 = LinearInterval("hg38", a.chr2, a.start2, a.end2)
        lin_ref_pairs.append((linRef1, linRef2))

    return lin_ref_pairs


def get_filtered_alignments(filename, offset=0):
    a = get_alignments(filename, offset)
    a = filter_alignments(a)
    return create_linear_intervals_from_alignments(a)

if __name__ == "__main__":
    a = get_alignments(DATA_PATH + "alignment.tmp")
    if DEBUG: print a
