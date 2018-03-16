#!/usr/bin/env python

'''
Convert the three-file raw Illumina data (forward, reverse, and index reads) into the
two-file format used by SmileTrain.
'''

import argparse, sys, os, itertools
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from Bio import SeqIO
import re

def new_id(location, direction):
    '''(location, barcode, direction) -> location/direction'''
    return "%s/%s" %(location, direction)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('good_align', help='good reads from mothur alignment')
    parser.add_argument('bed_in', help='bed files (tab)')
    parser.add_argument('one_out', help='output one reads fastq')
    
    args = parser.parse_args()
    #start an empty lookup list
    lookuplist = []
    printedlist = []
    
    f = open(args.bed_in, 'r')
    for line in f:
        new = re.split('\t', str(line))
        if new[3].endswith("/1"):
            new[3]  = new[3].rstrip("1")
            new[3] = new[3].rstrip("/")
        if new[3].endswith("/2"):
            new[3] = new[3].rstrip("2")
            new[3] = new[3].rstrip("/")
        lookuplist.append(new[3])

    f.close()

    with open(args.one_out, 'w') as go:
        for record in SeqIO.parse(args.good_align, "fasta"):
            new = re.split('_', str(record.id))
            beg = '_'.join(map(str, new[0:4]))
            end = ":".join(map(str, new[4:]))
            readid = ":".join(map(str, [beg, end]))
            # check that the reads match
            #if the ids are already formatted with /123, remove it
            if readid.endswith("/1"):
                readid  = readid.rstrip("1")
                readid = readid.rstrip("/")
            if readid.endswith("/2"):
                readid = readid.rstrip("2")
                readid = readid.rstrip("/")

            if readid in lookuplist:
                if readid not in printedlist:
                    # write the entries to their output files
                    printedlist.append(readid)
                    SeqIO.write(record, go, 'fasta')
 
