def divide_gen_file(genes_fn):
    # Divide into one file per chromosome
    f = open(genes_fn, "r")
    header = f.readline()
    from collections import defaultdict
    lines = defaultdict(list)
    for l in f.readlines():
        chrom = l.split()[2]
        if "_" in chrom:
            chrom = chrom.split("_")[0]
        lines[chrom].append(l)
    f.close()

    for chrom in lines:
        f = open("data/genes/genes_gencode_%s.txt" % chrom, "w")
        f.writelines([header])
        f.writelines(lines[chrom])
        f.close()

divide_gen_file("data/genes/genes_gencode.txt")