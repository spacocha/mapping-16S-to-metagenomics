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
    parser.add_argument('trans_in', help='file of which contigs are in the bins')
    parser.add_argument('log_in', help='which reads associated to which sequences')
    parser.add_argument('blast_in', help='blast list file')
    parser.add_argument('log_out', help='output table')
    parser.add_argument('fa_out', help='output fasta')
    
    args = parser.parse_args()
    bartrans={}
    lo=open(args.log_out, 'w')
    mi=open(args.meta_in, 'r')
    for line in mi:
        line=line.strip()
        new = re.split(',', str(line))
        bars= re.split('-', new[-1])
        if len(bars) == 2:
            newbar="%s~%s" %(bars[1], bars[0])
            bartrans[newbar]=new[1]
    mi.close()    

    transdict = {}
    ti=open(args.trans_in, 'r')
    for line in ti:
        line=line.strip()
        new = re.split('\t', str(line))
        transdict[new[1]] = new[0]

    ti.close()
    seqdict = {}
    fi=open(args.log_in, 'r')
    for line in fi:
        #this includes the base name, mothur name and alignment
        line=line.strip()
        new = re.split('\t', str(line))
        read1="%s/1" %(new[0])
        read2="%s/2" %(new[0])
        seqdict[read1]=new[2]
        seqdict[read2]=new[2]


    #now count how many blast hits associated with each bin
    countdict={}
    #also output a fasta file with the fastas for each bin
    #open blast list
    bi=open(args.blast_in, 'r')
    for line in bi:
        #remove \n
        line=line.strip()
        #split to get the barcode translation
        new = re.split('_', str(line))
        #only use it if there is a translations
        if new[2] in bartrans:
            #if so, open the actual report
            ri=open(line, 'r')
            for rline in ri:
                #remove any \n
                rline=rline.strip()
                #split into query, match, percent etc
                br = re.split('\t', str(rline))
                #only deal with reports that are in bins and hit with 100% id (blast has a prefilter of 1E-50)
                if br[1] in transdict:
                    if float(br[2]) == 100:
                        if br[0] in seqdict:
                            #these are the records that tell us how many hits per bin and which 16S sequence is associated with it
                            #contig match print(br[1])
                            #bin print(transdict[br[1]])
                            #sample print(bartrans[new[2]])
                            #print("%s\t%s\t%s\t%s" %(br[0], br[1], transdict[br[1]], seqdict[br[0]]))
                            #make the tuple to use as key
                            if transdict[br[1]] in countdict:
                                countdict[transdict[br[1]]]+=1
                            else:
                                countdict[transdict[br[1]]]=1
                            #print out the fasta sequence to a file for each bin
                            ioname="%s.%s.fa" %(args.fa_out, transdict[br[1]])
                            bo=open(ioname, 'a')
                            frecord = SeqRecord(Seq(seqdict[br[0]]), id=br[0])
                            SeqIO.write(frecord, bo, 'fasta')
                            bo.close()
    for keys in countdict:
	print("%s\t%s" %(keys, countdict[keys]))
	
