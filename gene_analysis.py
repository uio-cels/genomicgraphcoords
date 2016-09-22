from __future__ import print_function
from __future__ import absolute_import
from DbWrapper import DbWrapper
from offsetbasedgraph import LinearInterval
from main import get_flanks


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
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    alt_partitions = [par[1] for par in partitions]
    lin_ref_dict = dict(zip(names, alt_partitions))

    real_alt_genes = [gi for gi in gene_intervals if
                      lin_ref_dict[gi.chromosome][1].intersects(gi)]

    # Find genes that are partly in aligned regions
    intersecting_genes = [gi for gi in real_alt_genes if
                          gi.intersects(lin_ref_dict[gi.chromosome][0]) or
                          gi.intersects(lin_ref_dict[gi.chromosome][2])]

    print(len(intersecting_genes))
    print(len(real_alt_genes))
    print(len(intersecting_genes)/float(len(real_alt_genes)))


def get_overlap_code(gene, loci_partition):
    full_loci = loci_partition[0] + loci_partition[1] + loci_partition[2]
    return (
        not full_loci.contains(gene),
        gene.intersects(loci_partition[0]) or
        gene.intersects(loci_partition[2]),
        gene.intersects(loci_partition[1]),
        gene.intersects(loci_partition[0]) and
        gene.intersects(loci_partition[2]),
        gene.contains(full_loci))


def code_more(a, b):
    for v1, v2 in zip(a, b):
        if (not v1) and v2:
            return False
    return True


def code_less(a, b):
    for v1, v2 in zip(a, b):
        if v1 and not v2:
            return False
    return True


def code_equal(a, b):
    return a == b


def find_main_genes_by_codes(gene_intervals, lin_ref_dict, code,
                             cmpr=code_equal):
    """
    Find genes that overlap with an alt_loci in the way specified by code
    """
    return [gi for gi in gene_intervals if
            any([cmpr(get_overlap_code(gi, partition), code)
                 for partition in lin_ref_dict[gi.chromosome]])]


def get_alt_gene_stats():
    db = DbWrapper()
    genes = db.get_alt_genes()
    gene_intervals = [interval_from_gene(gene) for gene in genes]
    alt_loci_infos = db.get_alt_loci_infos(False)
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    alt_partitions = [par[1] for par in partitions]

    lin_ref_dict = {}
    for gene in gene_intervals:
        lin_ref_dict[gene.chromosome] = []

    for i, ali in enumerate(alt_loci_infos):
        chrom = ali["name"]
        if chrom not in lin_ref_dict:
            lin_ref_dict[chrom] = []
        lin_ref_dict[chrom].append(alt_partitions[i])

    pure_alt_code = (False, False, True, False, False)
    pure_flank_code = (False, True, False, False, False)
    alt_genes = find_main_genes_by_codes(
        gene_intervals, lin_ref_dict, pure_alt_code)
    flank_genes = find_main_genes_by_codes(
        gene_intervals, lin_ref_dict, pure_flank_code)
    spanning_genes = find_main_genes_by_codes(
        gene_intervals, lin_ref_dict, pure_alt_code, code_more)

    print("All genes", len(gene_intervals))
    print("Alt genes", len(alt_genes))
    print("Flank genes", len(flank_genes))
    print("Spanning genes", len(spanning_genes)-len(alt_genes))


def get_main_gene_stats():
    db = DbWrapper()
    genes = db.get_main_genes()
    gene_intervals = [interval_from_gene(gene) for gene in genes]
    alt_loci_infos = db.get_alt_loci_infos(False)
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    main_partitions = [par[0] for par in partitions]

    lin_ref_dict = {}
    for gene in gene_intervals:
        lin_ref_dict[gene.chromosome] = []

    for i, ali in enumerate(alt_loci_infos):
        chrom = ali["chrom"]
        if chrom not in lin_ref_dict:
            lin_ref_dict[chrom] = []
        lin_ref_dict[chrom].append(main_partitions[i])

    old_spanning_codes = (True, True, False, False, False)

    spanning_genes = find_main_genes_by_codes(
        gene_intervals, lin_ref_dict, old_spanning_codes)

    alt_code = (False, False, True, False, False)
    alt_genes = find_main_genes_by_codes(gene_intervals, lin_ref_dict,
                                         alt_code)

    new_spanning_genes = find_main_genes_by_codes(
        gene_intervals, lin_ref_dict, alt_code, code_more)

    pure_flankning_genes = find_main_genes_by_codes(
        gene_intervals, lin_ref_dict, (False, True, False, False, False))
    print("All genes", len(gene_intervals))
    print("Pure alt_genes", len(alt_genes))
    print("Spanning alt genes", len(new_spanning_genes)-len(alt_genes))
    print("Flanking genes", len(pure_flankning_genes))
    print("Flank spanning genes",  len(spanning_genes))


def calculate_original_main_spans():
    """
    Find all genes on the main path that intersects both
    with a flanking region and a pure main region.
    """
    db = DbWrapper()
    genes = db.get_main_genes()
    gene_intervals = [interval_from_gene(gene) for gene in genes]
    alt_loci_infos = db.get_alt_loci_infos(False)
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    main_partitions = [par[0] for par in partitions]

    lin_ref_dict = {}
    for gene in gene_intervals:
        lin_ref_dict[gene.chromosome] = []

    for i, ali in enumerate(alt_loci_infos):
        chrom = ali["chrom"]
        if chrom not in lin_ref_dict:
            lin_ref_dict[chrom] = []
        lin_ref_dict[chrom].append(main_partitions[i])

    intersecting_genes = [gi for gi in gene_intervals if
                          any([(gi.intersects(interval[0]) or
                                gi.intersects(interval[2])) and
                               not gi.intersects(interval[2]) for
                               interval in lin_ref_dict[gi.chromosome]])
                          ]
    crossing_genes = [gi for gi in intersecting_genes if
                      any([(interval[0].intersects(gi) and not interval[0].contains(gi)) or
                           (interval[2].intersects(gi) and not interval[2].contains(gi)) for
                           interval in lin_ref_dict[gi.chromosome]])]

    print(len(gene_intervals))
    print(len(intersecting_genes))
    print(len(crossing_genes))
    for gene in crossing_genes:
        print(gene.gene_name, gene)


def calculate_main_spans():
    """
    Find all genes on main paths that overlaps with real
    alt sequence, and find ratio of those that intersects both
    main and real alt.
    """
    db = DbWrapper()
    genes = db.get_main_genes()
    gene_intervals = [interval_from_gene(gene) for gene in genes]
    alt_loci_infos = db.get_alt_loci_infos(False)
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    main_partitions = [par[0] for par in partitions]
    lin_ref_dict = {}
    for gene in gene_intervals:
        lin_ref_dict[gene.chromosome] = []

    for i, ali in enumerate(alt_loci_infos):
        chrom = ali["chrom"]
        if chrom not in lin_ref_dict:
            lin_ref_dict[chrom] = []
        lin_ref_dict[chrom].append(main_partitions[i])

    print(len(genes), len(gene_intervals))

    real_alt_genes = [gi for gi in gene_intervals if
                      any([interval[1].intersects(gi) for
                           interval in lin_ref_dict[gi.chromosome]])]

    # Find genes that are partly in aligned regions
    intersecting_genes = [gi for gi in real_alt_genes if
                          any([gi.intersects(interval[0]) +
                               gi.intersects(interval[1]) +
                               gi.intersects(interval[2]) > 1 for
                               interval in lin_ref_dict[gi.chromosome]])
                          ]

    print(len(intersecting_genes))
    print(len(real_alt_genes))
    print(len(intersecting_genes)/float(len(real_alt_genes)))
    print(len(intersecting_genes)/float(len(genes)))


if __name__ == "__main__":
    get_alt_gene_stats()
    # get_main_gene_stats()

    # get_flanking_lins()
    # calculate_main_spans()
    exit(0)
