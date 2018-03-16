#!/usr/bin/env python

'''
Convert the three-file raw Illumina data (forward, reverse, and index reads) into the
two-file format used by SmileTrain.
'''

import argparse, sys, os, itertools
sys.path.append(os.path.normpath(os.path.abspath(__file__) + '/../..'))
from Bio import SeqIO
import re
import numpy

def new_id(location, direction):
    '''(location, barcode, direction) -> location/direction'''
    return "%s/%s" %(location, direction)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dist_in', help='distance files from mothur')
    parser.add_argument('one_out', help='output average and standard deviation')
    
    args = parser.parse_args()
    #start an empty lookup list

    lookuplist = []
    f = open(args.dist_in, 'r')
    o = open(args.one_out, 'w')
    for line in f:
        new = re.split(' ', str(line))
        new[-1] = new[-1].strip()
        lookuplist.append(float(new[2]))

    f.close()
    if len(lookuplist)> 1:
        ave= numpy.average(lookuplist, axis=0)
        std = numpy.std(lookuplist, axis=0)
        pstr="%f\t%f\t%d\n" %(ave, std, len(lookuplist))
        o.writelines(pstr)

 
