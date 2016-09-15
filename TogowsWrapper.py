"""
Simple methods for accessing the Togows API (http://togows.org/)
"""
import urllib2
from config import *
DEBUG = True

def get_sequence(loci_id, start=1, end=0):
    """
    Gets the sequence (as fasta format)
    of an alternative loci (e.g. chr1 or chr1_GL383519v1_alt)
    """
    url = "http://togows.org/api/ucsc/hg38/%s:%d-%d.fasta" % (loci_id, start, end)
    # if DEBUG: print "Fetching sequence from " + url
    try:
        sequence = urllib2.urlopen(url).read()
    except Exception as e:
        print url
        raise e
    return sequence


def save_sequence_to_fasta(loci_id, start, end, file_name = ""):
    """
    Collects a sequnce from the togows API, and saves it to file (fasta format)
    :param loci_id: The id of the loci (e.g. chr1 or chr11_KI270827v1_alt)
    :param start: Start position >= 1
    :param end: End position (inclusive (?))
    :param file_name: Alternative file name. If empty, one will be generated
    :return: Returns the file name
    """

    start += 1
    end += 1

    # Save url so that it can be presented later
    url = "http://togows.org/api/ucsc/hg38/%s:%d-%d.fasta" % (loci_id, start, end)
    import globals
    if "alt" in loci_id:
        globals.togows_alt_url = url
    else:
        globals.togows_main_url = url

    if file_name == "":
        file_name = DATA_PATH + "%s:%d-%d.fasta" % (loci_id, start, end)


    if DEBUG: print "=============== FETCHING SEQUENCE =============="
    import os.path
    if os.path.isfile(file_name):
        if DEBUG: print "File %s is chaced" % file_name
        return file_name

#    if DEBUG: print "###########"
#    if DEBUG: print file_name
#    if DEBUG: print loci_id, start, end
#    if DEBUG: print "---------------"

    curpath = os.path.abspath(os.curdir)
    if DEBUG: print curpath
    seq = get_sequence(loci_id, start, end)
    f = open(file_name, "w")
    f.write(seq)
    return file_name
