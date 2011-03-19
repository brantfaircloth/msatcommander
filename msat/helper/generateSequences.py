#!/usr/bin/env python
# encoding: utf-8
"""
generateSequences.py

Created by Brant Faircloth on 2009-03-28.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
import motif
import numpy

def generateSequence(motifs, outfile, max_rep = 14):
    for m in motifs:
        motif_length = m[0]
        for seq in m[1:]:
            prime5_seq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' #len=36
            prime3_seq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' #len=36
            #pdb.set_trace()
            rep_number = numpy.random.randint(6,max_rep)
            prime5_end = numpy.random.randint(16,36)
            prime3_end = numpy.random.randint(16,36)
            insert = seq*rep_number
            start = len(prime5_seq[:prime5_end])
            end = start + len(insert)
            seq_id = '>%s_%s_%s_(%s,%s)' % (motif_length, seq, rep_number, start, end)
            test_seq = prime5_seq[:prime5_end] + insert + prime3_seq[:prime3_end]
            outfile.write('%s\n%s\n' % (seq_id, test_seq))
            


def main():
    motifs = (motif.mononucleotide, motif.dinucleotide, motif.trinucleotide, motif.tetranucleotide, motif.pentanucleotide, motif.hexanucleotide)
    outfile = open('testRepeats.fa', 'w')
    generateSequence(motifs, outfile)
    outfile.close()


if __name__ == '__main__':
    main()

