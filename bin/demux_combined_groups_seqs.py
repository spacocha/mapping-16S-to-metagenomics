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
    parser.add_argument('fasta_in', help='combined_fasta file')
    parser.add_argument('groups_in', help='groups file')
    parser.add_argument('fasta_out', help='combined demux fasta file, like q.fst')
    groupsdict={}
    args = parser.parse_args()
    fi=open(args.groups_in, 'r')
    for line in fi:
        line=line.strip()
        new=re.split('\t', line)
        header=new[0]
        lib=new[1]
        groupsdict[str(header)]=str(lib)
    fi.close()

    count=1
    fo=open(args.fasta_out, 'w')
    for record in SeqIO.parse(args.fasta_in, "fasta"):
        if str(record.id) in groupsdict:
            record.id="sample=%s;%d/1" % (str(groupsdict[str(record.id)]), int(count))
            count+=1
            record.description=""
            SeqIO.write(record, fo, 'fasta')
        
