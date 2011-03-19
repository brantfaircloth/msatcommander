#!/usr/bin/env python
# encoding: utf-8
"""
msat/seqsearch.py

Part of 454_msatcommander.

Created by Brant Faircloth on 2009-03-29.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import re
import pdb
import string
from Bio.Seq import Seq

'''
Example use of class

>>> import seqsearch
>>> import motif
>>> from Bio.Seq import Seq

>>> dinucleotide = seqsearch.MicrosatelliteMotif(motif.dinucleotide, 6, \
False, name='Dinucleotide')
>>> seq = Seq('ACCTCGCGGGTGTCGCATATATATATATATATATATCATGTCTGCTCGGTGCGAGTATC')
>>> search = seqsearch.RegionSearch(seq)
>>> search.microsatellite(dinucleotide)
>>> search.matches
'''

class MicrosatelliteMotif:
    '''The class for all Motif search objects.
    Motif is initialized with:
    
        >>>msat = Motif(motif, repeat_length)
        
        :Parameters:
            - motif : a tuple containing the repeat motifs in which we are \
interested of the form (2, '(AC)', '(AG)', '(AT)', '(CG)') where the 0th \
element gives the motif length and the > 0th elements represent the actual\
motifs
            - min_length : the minimum length of a microsatellite of \
motif[element] that we are willing to accept
    '''
    def __init__(self, motif, min_length, perfect=False, name=None):
        '''Initialize a Microsatellite search instance
        
        :Parameters:
            - motif : a tuple containing the repeat motifs in which we are \
interested of the form (2, '(AC)', '(AG)', '(AT)', '(CG)') where the 0th \
element gives the motif length and the >0th elements represent the actual \
motifs
            - min_length : the minimum length of a microsatellite of \
motif[element] that we are willing to accept
            - perfect : optional flag ensuring only perfect repeats (no \
ambiguous bases) are found.  Default is False which ensures a repeat like \
ACACACNNACACAC is treated as a valid microsatellite.
            - name : optional value indicating the name of the search class
        '''
        self.motif = motif[1:]
        self.name = name
        self.repeat_length = motif[0]
        self.perfect = perfect
        if not self.perfect:
            # subtracting 1 so that we get >= min_length due to N*(?:%s)+ in 
            # regex
            self.min_length = min_length - 1
        else:
            # for perfect repeats, we don't have the wildcard
            self.min_length = min_length
        self._compile()
        
    def __repr__(self):
        if not self.name:
            return (('< MicrosatelliteMotif instance, min_length %s, \
repeat_length %s >') % (self.min_length, \
            self.repeat_length))
        else:
            return (('< MicrosatelliteMotif instance, min_length %s, \
repeat_length %s, name: %s >') % (self.min_length, \
            self.repeat_length, self.name))
    
    def _compile(self):
        '''Compiled regular expressions for the given motif (PRIVATE)
        
        Called by __init__ on instantiation.This is primarily to ensure (1) \
each set of regexes is generated 1X, (2) we are not wasting time generating \
regexes for unnecessary repeat motifs and (3) we are creating persistent \
search classes'''
        self.compiled = ()
        for f in self.motif:
            bases = string.maketrans('AGCTagct','TCGAtcga')
            r = string.translate(f, bases)[::-1]
            if not self.perfect:
                # forward motif wild-card
                wild1 = 'N*(?:%s)+' % f
                # reverse complement wild-card
                wild2 = 'N*(?:%s)+' % r
                rep_match = '(?:%s){%s,}%s|(?:%s){%s,}%s' % (f, \
                self.min_length, wild1, r, self.min_length, wild2)
            else:
                rep_match = '(?:%s){%s,}|(?:%s){%s,}' % (f, self.min_length, \
                r, self.min_length)
            rx = re.compile(rep_match,re.IGNORECASE)
            self.compiled += (rx,)
        
class RegionSearch:
    '''Base class for all RegionSearch objects.
    RegionSearch is initialized with:
    
        >>> search = RegionSearch(bioSeq)
        
        :Parameters:
        
            -bioSeq : a Biopython Seq instance
    '''
    def __init__(self, bioSeq):
        
        assert isinstance(bioSeq,Seq), 'This is not a Biopython Seq instance'
        self.forward = str(bioSeq)
        self.length = len(self.forward)
        self.matches = {}
        
    def __repr__(self):
        '''Keep representation similar to __repr__ of Bio.Seq

        Code borrowed and modified from Bio.Seq.Seq
        '''

        if len(self.forward) > 60 :
            return "%s('%s...%s')" % (self.__class__.__name__, \
            self.forward[:54], self.forward[-3:])
        else :
            return "%s(%s)" % (self.__class__.__name__, self.forward)
            
    def _generalized_search(self, seq, msat):
        '''Generalized microsatellite search function (PRIVATE)
        
        Receives input from self.microsatellite
        '''
        for repeat in range(len(msat.compiled)):
            temp_match = ()
            for m in msat.compiled[repeat].finditer(seq):
                temp_match += ((m.span(),m.span()[0],self.length-m.span()[1]),)
            if temp_match:
                self.matches[msat.motif[repeat]] = temp_match
                
    def microsatellite(self, msat):
        '''Sequence instance search method
        
        :Parameters:
        
            -msat : a Motif instance
        
        We are now only searching the upper sequence string (5'->3'), but \
using both the forward and reverse complement of the lowest alphabetical, \
non-complementary microsat source.  This rids us of length errors that you \
tend to see with certain microsat motifs differing btw. the upper and lower \
strands (i,e. non-palindromic & overlapping).
        '''
        self._generalized_search(self.forward, msat)