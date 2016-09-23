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


def create_lin_ref_dict(gene_intervals, alt_loci_infos,
                        partitions, key="chrom"):
    lin_ref_dict = {}
    for gene in gene_intervals:
        lin_ref_dict[gene.chromosome] = []

    for i, ali in enumerate(alt_loci_infos):
        chrom = ali[key]
        if chrom not in lin_ref_dict:
            lin_ref_dict[chrom] = []
        lin_ref_dict[chrom].append(partitions[i])

    return lin_ref_dict


def is_gene_equal(main_gene, alt_gene, partition):
    main_partition = partition[0]
    alt_partition = partition[1]
    if main_partition[0].intersects(main_gene):
        # print(main_gene, "first")
        main_start = main_gene.start-main_partition[0].start
        main_end = main_gene.end-main_partition[0].start
        # print(main_start, alt_gene.start)
        # print(main_end, alt_gene.end)
        t = main_start == alt_gene.start
        return t and main_end == alt_gene.end

    if main_partition[2].intersects(main_gene):
        # print(main_gene, "last")
        main_start = main_partition[2].end-main_gene.start
        alt_start = alt_partition[2].end-alt_gene.start
        main_end = main_partition[2].end-main_gene.end
        alt_end = alt_partition[2].end-alt_gene.end
        # print(main_start, alt_start)
        # print(main_end, alt_end)
        t = main_start == alt_start
        return t and main_end == alt_end

    return False


def check_flanking_equality(main_genes, alt_genes, alt_loci_infos):
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    lin_ref_dict_alt = create_lin_ref_dict(
        alt_genes, alt_loci_infos, [p[1] for p in partitions], "name")
    lin_ref_dict_main = create_lin_ref_dict(
        main_genes, alt_loci_infos, [p[0] for p in partitions])
    flanking_code = (False, True, False, False, False)

    main_flank_genes = find_main_genes_by_codes(main_genes, lin_ref_dict_main, flanking_code)
    alt_flank_genes = find_main_genes_by_codes(alt_genes, lin_ref_dict_alt, flanking_code)
    partition_map = {ali["name"]: partition for (ali, partition) in
                     zip(alt_loci_infos, partitions)}

    gene_pairs = []
    solo_genes = []
    for alt_flank_gene in alt_flank_genes:
        partition = partition_map[alt_flank_gene.chromosome]
        for gene in main_flank_genes:
            if is_gene_equal(gene, alt_flank_gene, partition):
                gene_pairs.append((alt_flank_gene, gene))
                # print("Found equal: (%s, %s)" % (alt_flank_gene, gene))
                break
        else:
            solo_genes.append(alt_flank_gene)

    main_gene_pairs = []
    main_solo_genes = []
    for main_flank_gene in main_flank_genes:
        partition = [p for p in partitions if
                     p[0][0].contains(main_flank_gene) or
                     p[0][2].contains(main_flank_gene)][0]
        for alt_gene in alt_flank_genes:
            if is_gene_equal(main_flank_gene, alt_gene, partition):
                main_gene_pairs.append((main_flank_gene, alt_gene))
                break
        else:
            main_solo_genes.append(main_flank_gene)
    print(len(main_gene_pairs))
    print(len(main_solo_genes))
    # print("%s has no match\n %s" % (alt_flank_gene, partition[0][0]))
    # print("\n".join([(p[0].gene_name + " " + p[1].gene_name)
    # for p in gene_pairs]))
    print("\n".join([s.gene_name for s in solo_genes]))
    print("-----------------------")
    print("\n".join([s.gene_name for s in main_solo_genes]))
    print(len(gene_pairs))
    print(len(solo_genes))


def get_main_gene_stats():
    db = DbWrapper()
    genes = db.get_main_genes()
    alt_genes = db.get_alt_genes()
    gene_intervals = [interval_from_gene(gene) for gene in genes]
    alt_gene_intervals = [interval_from_gene(gene) for gene in alt_genes]
    alt_loci_infos = db.get_alt_loci_infos(False)
    return check_flanking_equality(gene_intervals, alt_gene_intervals, alt_loci_infos)
    
    partitions = [get_flanks(ali) for ali in alt_loci_infos]
    main_partitions = [par[0] for par in partitions]
    alt_partitions = [par[1] for par in partitions]

    lin_ref_dict = {}
    for gene in gene_intervals:
        lin_ref_dict[gene.chromosome] = []

    partition_dict = {}

    for i, ali in enumerate(alt_loci_infos):
        chrom = ali["chrom"]
        if chrom not in lin_ref_dict:
            lin_ref_dict[chrom] = []
        lin_ref_dict[chrom].append(main_partitions[i])
        partition_dict[main_partitions[i]] = alt_partitions[i]

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
    # get_alt_gene_stats()
    get_main_gene_stats()

    # get_flanking_lins()
    # calculate_main_spans()
    exit(0)
