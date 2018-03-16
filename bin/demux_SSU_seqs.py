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
import subprocess
from subprocess import Popen, PIPE, STDOUT

def new_id(location, direction):
    '''(location, barcode, direction) -> location/direction'''
    return "%s/%s" %(location, direction)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_in', help='list of files')
    parser.add_argument('fasta_out', help='combined demux fasta file, like q.fst')
    
    args = parser.parse_args()
    fi=open(args.fasta_in, 'r')
    fo=open(args.fasta_out, 'w')
    count=1
    for line in fi:
        line=line.strip()
        new=re.search('prefix.+.ssu.out', line)
        print(new.group(0))
        for record in SeqIO.parse(line, "fasta"):
            record.id="sample=%s;%d/1" % (str(new.group(0)), int(count))
            count+=1
            record.description=""
            SeqIO.write(record, fo, 'fasta')
