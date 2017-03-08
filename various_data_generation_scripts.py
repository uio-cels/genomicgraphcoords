"""
Contains various methods for generating/processing data.
"""

sizes = {}
f = open("grch38.chrom.sizes")
for line in f.readlines():
    l = line.split()
    sizes[l[0]] = int(l[1])


def create_alt_loci_file():
    f = open("grch38_alternative_loci.txt")
    lines_out = []
    for line in f.readlines():
        if line.startswith("#"):
            continue

        l = line.split()
        if l[4] != "alt-scaffold":
            continue

        chr = l[1]
        start = l[2]
        stop = l[3]
        id = "chr" + chr + "_" + l[5].replace(".", "v") + "_alt"
        size = sizes[id]

        lines_out.append("%s    chr%s  %s  %s %d\n" % (id, chr, start, stop, size))

    f.close()
    f2 = open("grch38_alt_loci.txt", "w")
    f2.writelines(lines_out)
    f2.close()


def create_alt_loci_file_from_db():
    f = open("hgTables.txt")
    lines_out = []
    for line in f.readlines():
        if line.startswith("#"):
            continue

        l = line.split()
        if "alt" in l[1]:
            continue

        chr = l[1]
        start = int(l[2]) + 1
        stop = l[3]
        id = l[4]
        size = sizes[id]

        lines_out.append("%s    %s  %s  %s %d\n" % (id, chr, start, stop, size))

    f.close()
    f2 = open("grch38_alt_loci_from_db.txt", "w")
    f2.writelines(lines_out)
    f2.close()

def get_genes(alt_loci, chr, start, end):
    # writes all genes on the alt loci and between start and end on the main chromosome to file
    fn = "genes_" + alt_loci + ".txt"
    out = []
    f = open("genes.txt")
    for line in f.readlines():
        d = line.split()
        if d[1] == chr and int(d[3]) >= start and int(d[4]) <= end:
            out.append(line)
        elif d[1] == alt_loci:
            out.append(line)

    o = open(fn, "w")
    o.writelines(out)
    o.close()
    f.close()

#get_genes("chr13_KI270838v1_alt", "chr13", 112182875, 112505943)


def curate_alignment_files():
    """
    Renames all alignment files downloaded from ncbi
    """

    # Find mapping from id to chrom
    alts = {}
    f = open("grch38.chrom.sizes")
    for line in f.readlines():
        l = line.split()
        id = l[0]
        if "alt" in id:
            s = id.split("_")
            chrom = s[0]
            alts[id.split("_")[1]] = chrom


    alignments_dir = "/home/ivar/alignments"
    import glob
    files = glob.glob("%s/*.gff" % alignments_dir)

    for fname in files:
        with open(fname, "r") as f:
            alt_id = fname.split("/")[-1].replace("gff", "").split("_")[0].replace(".", "v")
            chrom = alts[alt_id]
            alt_id = chrom + "_" + alt_id
            #print(alt_id)

            # Find first non-comment line
            cigar_line = ""
            for l in f.readlines():
                if not l.startswith("#"):
                    cigar_line = l
                    break

            #print(cigar_line)
            assert cigar_line != ""

            f2 = open(alignments_dir + "/" + alt_id + "_alt.alignment", "w")
            text = cigar_line
            #print(fname)
            s = text.split("Gap=")
            if len(s) < 2:
                print("No alignments for %s " % fname)
                continue
            cigar = s[1].split("#")[0]

            # Find main chr start/end
            s = cigar_line.split()
            #print(s)
            main_start = s[3]
            main_end = s[4]
            # Find alt pos
            s = text.split("Target=")[1].split(";")[0].split()
            alt_start = s[1]
            alt_end = s[2]

            line = "%s,%s,%s,%s,%s" % (main_start, main_end, alt_start, alt_end, cigar)
            print("%s,%s,%s,%s,%s" % (alt_id, main_start, main_end, alt_start, alt_end))
            f2.write(line)

            #print(cigar)
            f2.close()

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
        f = open("genes/genes_refseq_%s.txt" % chrom, "w")
        f.writelines([header])
        f.writelines(lines[chrom])
        f.close()

divide_gen_file("genes_refseq.txt")

#curate_alignment_files()

# create_alt_loci_file_from_db()