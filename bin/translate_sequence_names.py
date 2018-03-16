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
    parser.add_argument('meta_in', help='metadata file to map file names to samples')
    parser.add_argument('fasta_in', help='fasta file list')
    parser.add_argument('prefix_out', help='output prefix')
    
    args = parser.parse_args()
    bartrans={}
    mi=open(args.meta_in, 'r')
    for line in mi:
        line=line.strip()
        new = re.split(',', str(line))
        bars= re.split('-', new[-1])
        if len(bars) == 2:
            newbar="%s~%s" %(bars[1], bars[0])
            bartrans[newbar]=new[1]
    mi.close()    

    fi=open(args.fasta_in, 'r')
    for line in fi:
        #remove \n
        line=line.strip()
        new = re.split('_', str(line))
        if new[2] in bartrans:
            outputfile="%s_%s.fa" % (args.prefix_out, bartrans[new[2]])
            count = SeqIO.convert(line, "fastq", outputfile , "fasta")
            print("Converted %i records" % count)
