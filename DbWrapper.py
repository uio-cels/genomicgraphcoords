from __future__ import print_function
from __future__ import absolute_import
from builtins import object
import pymysql
from config import *
import numpy as np
import pickle

class Gene(object):
    def __init__(self, chrom_id, start, end, exon_starts, exon_ends):
        self.chrom_id = chrom_id
        self.start = start
        self.end = end
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends


class DbWrapper(object):
    """
    Abstraction layer for collecting genomic data from the UCSC database.
    """

    def __init__(self):
        self.chrom_lengths = {}
        self.db_ucsc = self.connect_ucsc()

        # Create a mapper from alt loci ID to new alt loci id
        self._create_loci_id_mapper()

    def _create_loci_id_mapper(self):
        self.alt_loci_names = {}
        mapping_file = DATA_PATH + "../hg38_alt_loci.bed"
        f = open(mapping_file)
        for l in f.readlines():
            if l[0] != "#":
                d = l.split()
                id = d[5]
                id = id.replace(".", "v")
                id += "_alt"
                id =  "chr" + d[1] + "_" + id
                self.alt_loci_names[id] = d[0] + "." + d[-1][-1:]

    def alt_loci_pretty_name(self, name):
        if name in self.alt_loci_names:
            return self.alt_loci_names[name]
        else:
            return name

    def connect_ucsc(self):
        return pymysql.connect(host="genome-mysql.cse.ucsc.edu",
                               database="hg38", user="genome")

    def fetch_all(self, query):
        """
        Returns all fetched rows from the given query as a list.
        :param query: The mysql query string
        :return: Returns a list of dicts (one list element per row)
        """

        # Check if query is cached in a text file. Saves UCSC for unecessary
        # databse calls.
        cached_file = DATA_PATH + "cached_query_" + self._md5(query)
        if os.path.isfile(cached_file):
            f = open(cached_file, 'rb')
            cached_result = pickle.load(f)
            f.close()
            return cached_result

        cursor = self.db_ucsc.cursor(pymysql.cursors.DictCursor)
        cursor.execute(query)
        rows = cursor.fetchall()
        cursor.close()

        # Cache the result to file
        afile = open(cached_file, 'wb')
        pickle.dump(rows, afile)
        afile.close()

        return rows

    def get_alt_loci_infos(self):
        """
        Gets all the alternative loci
        :return: A list of dicts (all alternative loci)
        """
        alt_loci = []

        res = self.fetch_all(
            r"SELECT alt.*, (select size from chromInfo i where i.chrom = alt.name) as length  FROM altLocations alt where name LIKE '%\_%'")
        if len(res) == 0:
            raise Exception("No results from database on alt loci")

        # Don return overlapping alt loci

        for alt1 in res:
            is_overlapping = False
            for alt2 in alt_loci:
                if alt1["chrom"] != alt2["chrom"]:
                    continue

                if alt1["chromStart"] < alt2["chromEnd"] and \
                            alt1["chromEnd"] > alt2["chromStart"] or \
                            alt2["chromStart"] < alt1["chromEnd"] and \
                            alt2["chromEnd"] > alt1["chromStart"]:
                    is_overlapping = True
                    break

            if not is_overlapping:
                alt_loci.append(alt1)

#        if DEBUG: print len(alt_loci)
#        if DEBUG: print len(res)
        return alt_loci

    def alt_loci_info(self, alt_loci_id):
        """
        Returns information about an alt loci, i.e the relative position
        on the main path in the genome ()
        :param alt_loci_id: the id of the alt loci, e.g. chr11_KI270827v1_alt
        :return: Returns a dict containing (chrom, chromStart, chromEnd)
        """
        res = self.fetch_all("SELECT chrom, chromStart, chromEnd, name FROM altLocations where name='%s'" % alt_loci_id)
        if len(res) == 0:
            raise Exception("No results from database on alt loci %s. This alt loci may not exist" % alt_loci_id)

        res2 = self.fetch_all("SELECT chromEnd as length from altLocations where chrom = '%s'" % alt_loci_id)
        res[0]["length"] = res2[0]["length"]
        return res[0]

    def get_chrom_lengths(self):
        # Gets all chromosome lengths
        res = self.fetch_all("SELECT chrom, size FROM chromInfo")
        for r in res:
            self.chrom_lengths[r["chrom"]] = r["size"]

    def chrom_length(self, chromosome_id):
        if chromosome_id in self.chrom_lengths:
            return self.chrom_lengths[chromosome_id]
        else:
            self.get_chrom_lengths()
            return self.chrom_lengths[chromosome_id]

    def alt_loci_names(self):
        res = self.fetch_all(
            r"SELECT chrom FROM altLocations where chrom LIKE '%\_%'")
        return [r["chrom"] for r in res]

    def genes_crossing_position(self, chrom_id, position):
        """
        Returns all genes starting before and ending after the given position
        :param chrom_id: alt locus or chrom id
        :param position: offset
        :return: Returns a list of dicts (each element representing one gene)
        """
        query = r"SELECT k.*, r.geneSymbol as gname FROM knownGene k, kgXref r where r.kgID = k.name AND k.chrom LIKE '%s' and k.txStart < %d and k.txEnd > %d" % (chrom_id, position, position)
        res = self.fetch_all(query)
        return res

    def get_gene(self, gene_id):
        query = r"SELECT * FROM knownGene where name='%s"

    def _md5(self, string):
        # Used to store cached queries, by hashing the query string
        import hashlib
        m = hashlib.md5()
        m.update(string.encode('utf-8'))
        return m.hexdigest()

if __name__ == "__main__":
    db = DbWrapper()
    #res = db.alt_loci_info("chr11_KI270827v1_alt")
    if DEBUG: print(db.genes_crossing_position("chr1", 100000))
