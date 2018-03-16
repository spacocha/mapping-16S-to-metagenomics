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
    parser.add_argument('fasta_in', help='fasta file of 16S sequences from SSUsearch')
    parser.add_argument('fastq_in', help='list of raw read files')
    parser.add_argument('mat_out', help='output OTU table')
    parser.add_argument('fa_out', help='output fasta file')
    
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
    count=1
    for record in SeqIO.parse(args.fasta_in, "fasta"):
        OTUid=str("seq%d" % count)
        print(OTUid)
        if record.seq not in fadict:
            newrecord = SeqRecord(Seq(str(record.seq)), id=OTUid)
            fadict[str(record.seq)]=OTUid
            count+=1
            print("Newrecord %s %s" %(str(record.seq), str(OTUid)))
            printlist.append(newrecord)
           
    SeqIO.write(printlist, args.fa_out, "fasta")

    OTUmatdict={}
    sampledict={}
    qi = open(args.fastq_in, 'r')
    for line in qi:
        line=line.strip()
        new = re.split('_', str(line))
        #check that we know which sample it came from
        print (new[2])
        if new[2] in bartrans:
            for seq in fadict:
                p=re.compile(str(seq))
                for record in SeqIO.parse(line, "fastq"):
                    #it could also be the reverse complement
                    if p.match(str(record.seq)):
                        #add one to the sample in bartrans
                        OTUmatdict[fadict[seq]][sample]+=1
                        sampledict[sample]+=1


    qi.close()

    with open(args.mat_out, 'w') as mo:
        #unpack the OTUmatdict headers
        mo.writelines("OTU\t")
        for sample in sampledict:
            mo.writelines(sample)
        mo.writelines("\n")

        #unpack OTUmatdict OTUs
        for seq in OTUmatdict:
            mo.writelines("%s" % seq)
            for sample in sampledict:
                if OTUmatdict[seq][sample]:
                    mo.writelines("\t%s" % OTUmatdict[seq][sample])
                else:
                    mo.writelines("\t0")
            mo.writelines("\n")

    mo.close()
