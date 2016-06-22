class LinearInterval(object):
    def __init__(self, genome_id, chromosome, start, end, strand="+"):
        self.genome_id = genome_id
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.gene_name = ""  # Used in some cases

    def length(self):
        return self.end-self.start

    def contains(self, other):
        if not (self.genome_id == other.genome_id):
            return False
        if not self.chromosome == other.chromosome:
            return False
        if self.start > other.start or self.end < other.end:
            return False
        return True

    def intersects(self, other):
        if not (self.genome_id == other.genome_id):
            return False
        if not self.chromosome == other.chromosome:
            return False
        if (self.start < other.end and self.end > other.start) or \
                (other.start < self.end and other.end > self.start):
            return True

        return False
        """
        if (self.start > other.start) == (self.end > other.end):
            return False
        return True
        """

    def __repr__(self):
        return "Lin seg in species %s, %s [%d, %d] %s" % \
               (self.genome_id, self.chromosome,
                self.start, self.end, self.strand)

    def label(self):
        return "%s" % (self.chromosome)

    def __str__(self):
        return "%s from %d % d on chr %s" % (self.label(), self.start, self.end, self.chromosome)
