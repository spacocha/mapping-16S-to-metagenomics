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
    parser.add_argument('one_in', help='input one reads fastq')
    parser.add_argument('two_in', help='input two reads fastq')
    parser.add_argument('one_out', help='output one reads fastq')
    parser.add_argument('two_out', help='output two reads fastq')
    parser.add_argument('log_out', help='output for sequence maps')
    
    args = parser.parse_args()
    #start an empyt lookup list
    list= []

    lo=open(args.log_out, 'w')
    for record in SeqIO.parse(args.good_align, "fasta"):
        new = re.split('_', str(record.id))
        beg = '_'.join(map(str, new[0:4]))
        end = ":".join(map(str, new[4:]))
        readid = ":".join(map(str, [beg, end]))
        if readid.endswith("/1"):
            readid = readid.rstrip("1")
            readid = readid.rstrip("/")
        if readid.endswith("/2"):
	    readid = readid.rstrip("2")
            readid = readid.rstrip("/")
        #Make this into a lookup table
        list.append(readid)
        trans=[readid, str(record.id), str(record.seq)] 
        pstr="\t".join(map(str, trans))
        lo.writelines(pstr)
        lo.writelines("\n")
    lo.close()    

    with open(args.one_out, 'w') as fo, open(args.two_out, 'w') as ro:
        fi = SeqIO.parse(args.one_in, 'fastq')
        si = SeqIO.parse(args.two_in, 'fastq')
        
        #First I need to make a lookup table of all the headers I want to keep
            
        for fr, sr in itertools.izip(fi, si):
            # check that the reads match
            #if the ids are already formatted with /123, remove it
            if fr.id.endswith("/1"):
                fr.id  = fr.id.rstrip("1")
                fr.id = fr.id.rstrip("/")
            if sr.id.endswith("/2"):
                sr.id = sr.id.rstrip("2")
                sr.id = sr.id.rstrip("/")
                
            if not (fr.id == sr.id):
                raise RuntimeError("ids in corresponding entries did not match: %s %s %s" %(fr.id, rr.id))
            
            # only keep 
            if fr.id in list:                
                # amend the ids on the forward and reverse entries
                fr.id = new_id(fr.id, '1')
                sr.id = new_id(sr.id, '2')
                fr.description = ''
                sr.description = ''
            
                # write the entries to their output files
                SeqIO.write(fr, fo, 'fastq')
                SeqIO.write(sr, ro, 'fastq')
