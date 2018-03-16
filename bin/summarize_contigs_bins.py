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
    parser.add_argument('one_in', help='list of fasta files to map')
    parser.add_argument('log_out', help='output for sequence maps')
    
    args = parser.parse_args()
    #start an empyt lookup list

    lo=open(args.log_out, 'w')
    io=open(args.one_in, 'r')
    for line in io:
        line = line.strip()
        for record in SeqIO.parse(line, "fasta"):
            trans=[line, str(record.id)] 
            pstr="\t".join(map(str, trans))
            lo.writelines(pstr)
            lo.writelines("\n")
    lo.close()    

