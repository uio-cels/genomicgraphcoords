from DbWrapper import DbWrapper
from LinearInterval import LinearInterval
from main import get_flanking_alignments


def interval_from_gene(gene):
    interval = LinearInterval("hg38",
                              gene["chrom"],
                              gene["txStart"],
                              gene["txEnd"])
    interval.gene_name = gene["name"]
    return interval


def get_flanking_lins():
    db = DbWrapper()
    genes = db.get_alt_genes()
    gene_intervals = map(interval_from_gene,
                         genes)

    alt_loci_infos = db.get_alt_loci_infos(False)
    names = map(lambda x: x["name"], alt_loci_infos)
    alignments = map(lambda x: get_flanking_alignments(x["name"], x),
                     alt_loci_infos)
    lin_refs = map(lambda x: [y[1] for y in x], alignments)
    lin_ref_dict = dict(zip(names, lin_refs))

    # Remove genes that are entierly inside aligned regions
    not_contains = filter(
        lambda gi: not any([y.contains(gi) for y in lin_ref_dict[gi.chromosome]]),
        gene_intervals)

    # Find genes that are partly in aligned regions

    intersects = filter(
        lambda gi: any([gi.intersects(y) for y in lin_ref_dict[gi.chromosome]]),
        not_contains)

    print len(intersects)
    print len(not_contains)
    print len(intersects)/float(len(not_contains))


if __name__ == "__main__":
    get_flanking_lins()
    exit(0)
