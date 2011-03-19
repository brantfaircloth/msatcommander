#!/usr/bin/env python
# encoding: utf-8
"""
msat/finder.py

Part of 454_msatcommander.

Created by Brant Faircloth on 2009-03-29.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""
import pdb
import motif
import seqsearch

def createMotifInstances(motif, min_length, perfect):
    return seqsearch.MicrosatelliteMotif(motif, min_length, perfect)

def genMotifCollection(options):
    possible_motifs = (motif.mononucleotide, motif.dinucleotide, \
    motif.trinucleotide, motif.tetranucleotide, motif.pentanucleotide, \
    motif.hexanucleotide)
    # add optional lengths
    possible_motifs = zip(options.min_length, possible_motifs)
    motif_collection = ()
    if options.scan_type == 'all':
        for m in possible_motifs:
            motif_collection += (createMotifInstances(m[1], m[0], \
            options.perfect),)
    elif '+' in options.scan_type:
        # subtracting 1 so that we get >= options.scan_type
        scan = int(options.scan_type[0]) - 1
        for m in possible_motifs[scan:]:
            motif_collection += (createMotifInstances(m[1], m[0], \
            options.perfect),)
    elif '-' in options.scan_type:
        scan_start = int(options.scan_type[0]) - 1
        scan_stop = int(options.scan_type[2])
        for m in possible_motifs[scan_start:scan_stop]:
            motif_collection += (createMotifInstances(m[1], m[0], \
            options.perfect),)
    else:
        # no iteration here because tuple != nested
        scan = int(options.scan_type[0]) - 1
        motif_collection += (createMotifInstances(possible_motifs[scan][1], \
        possible_motifs[scan][0], options.perfect),)
    return motif_collection

def stdOut(name, matches):
    if matches:
        for msat in matches:
            for match in matches[msat]:
                start = match[0]
                end = match[1]
                length = (end - start) / len(msat)
                print '%s:\t(%s)^%s\trepeat between bases\t%s\tand\t%s' % \
                (name, msat, length, start, end)
    else:
        print '%s:\t No microsatellites found (or repeat region < \
min_length)' % name