#!/usr/bin/env python

'''
Convert the three-file raw Illumina data (forward, reverse, and index reads) into the
two-file format used by SmileTrain.
'''

import argparse, sys, os, itertools
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re

def new_id(location, direction):
    '''(location, barcode, direction) -> location/direction'''
    return "%s/%s" %(location, direction)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_in', help='fasta file of 16S sequences from SSUsearch')
    parser.add_argument('fa_out', help='output fasta file')
    
    args = parser.parse_args()

    fadict={}
    printlist=[]
    count=1
    for record in SeqIO.parse(args.fasta_in, "fasta"):
        if record.seq not in fadict:
            OTUid=str("seq%d" % count)
            newrecord = SeqRecord(Seq(str(record.seq)), id=OTUid)
            fadict[str(record.seq)]=OTUid
            count+=1
            printlist.append(newrecord)
           
    SeqIO.write(printlist, args.fa_out, "fasta")

