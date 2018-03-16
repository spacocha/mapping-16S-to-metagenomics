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
    parser.add_argument('meta_in', help='metadata file to map file names to samples')
    parser.add_argument('fasta_in', help='dereplicated fasta')
    parser.add_argument('fastq_in', help='list of raw read files')
    parser.add_argument('mat_out', help='output OTU table')
    
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

    fadict={}
    printlist=[]
    cmd=[]
    count=1
    for record in SeqIO.parse(args.fasta_in, "fasta"):
            fadict[str(record.seq)]=str(record.id)
           
    OTUlist=[]
    samplelist=[]
    seqlist=[]
    qi = open(args.fastq_in, 'r')
    for line in qi:
        line=line.strip()
        new = re.split('_', str(line))
        #check that we know which sample it came from
        if new[2] in bartrans:
            #clear sample list
            samplelist.append(bartrans[new[2]])
            for seq in fadict:
                if fadict[seq] not in seqlist:
                    #add this to the list of sequences to process
                    seqlist.append(fadict[seq])
                #system call and use grep to get the total number of times that the seuqence is found in the file
                cmd="grep -c %s %s" % (seq, str(line))
                p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
                output = p.stdout.read()
                samplelist.append(bartrans[new[2]])
                print "%s\t%s\t%s\t%d" % (fadict[seq], bartrans[new[2]], str(line), int(output))
                #put this into some matrix format to recall later
            

    qi.close()

