#!/usr/bin/env python
# encoding: utf-8
"""
seqsearch_tests.py

Created by Brant Faircloth on 2009-03-27.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import seqsearch
import pdb
import motif
from nose import tools
from Bio.Seq import Seq
from Bio import SeqIO

class Test_1_Motif:
    '''Functions for testing our motifs'''
    def test_1_nose(self):
        '''[Sanity] : ensure our testing environment is fundamentally sane'''
        assert 4 == 4
        assert 'a' == 'a'
        assert 5678 - 2 == 5676
    def test_2_mono(self):
        '''[Motif] : ensure that mononucleotide repeats are sane'''
        assert motif.mononucleotide[1:] == ('A','C')
    def test_3_dinuc(self):
        '''[Motif] : ensure that dinucleotide repeats are sane'''
        assert motif.dinucleotide[1:] == ('AC','AG','AT','CG')
    def test_4_trinuc(self):
        '''[Motif] : ensure that trinucleotide repeats are sane'''
        assert motif.trinucleotide[1:] == ('AAC','AAG','AAT','ACC','ACG','ACT','AGC','AGG','ATC','CCG')
    def test_5_trinuc_ne(self):
        '''[Motif] : ensure that trinucleotide repeats don't match what they're not supposed to'''
        assert motif.trinucleotide is not ('AAAC')

class Test_3_CompiledMotifs:
    '''Functions for testing to ensure motifs are compiled'''
    def test_1_AC_repeat(self):
        '''[MicrosatelliteMotif] : ensure 1 perfect repeat area is found and that the length is correct'''
        # base seq (AC)12
        base_seq = 'ACAAACAGAGAAATATGACACACACACACACACACACACACGGTGTTGATAGTCAGGTGA'
        dinuc = seqsearch.MicrosatelliteMotif(motif.dinucleotide, 6, True, 'Dinucleotide')
        ac_repeat = dinuc.compiled[0]
        matches = ac_repeat.finditer(base_seq)
        i = 0
        for m in matches:
            b = m.span()
            assert (b[1]-b[0])/2 == 12
            i += 1
        assert i == 1
        
    def test_2_AC_repeat(self):
        '''[MicrosatelliteMotif] : ensure 1 imperfect repeat area is found and that the length is correct'''
        # base seq (AC)12
        base_seq = 'ACAAACAGAGAAATATGACACACACACACACACACACACNNACGGTGTTGATAGTCAGGTGA'
        dinuc = seqsearch.MicrosatelliteMotif(motif.dinucleotide, 6, False, 'Dinucleotide')
        ac_repeat = dinuc.compiled[0]
        matches = ac_repeat.finditer(base_seq)
        i = 0
        for m in matches:
            b = m.span()
            assert (b[1]-b[0])/2 == 13
            i += 1
        assert i == 1
    
    def test_3_AC_multiple(self):
        '''[MicrosatelliteMotif] : ensure that multiple perfect repeats are found and both are of length'''
        base_seq = 'ACAAACAGAGAAATATGACACACACACACACACACACACACGGTGTTGATAGTCAGGTGAACAAACAGAGAAATATGACACACACACACACACACACACACGGTGTTGATAGTCAGGTGA'
        dinuc = seqsearch.MicrosatelliteMotif(motif.dinucleotide, 6, True, 'Dinucleotide')
        ac_repeat = dinuc.compiled[0]
        matches = ac_repeat.finditer(base_seq)
        i = 0
        #tools.set_trace()
        for m in matches:
            b = m.span()
            print b
            assert (b[1]-b[0])/2 == 12, 'Repeat lengths are not as expected'
            i += 1
        assert i == 2
        
    def test_4_AC_short(self):
        '''[MicrosatelliteMotif] : ensure that short, perfect repeats ARE NOT found'''
        base_seq = 'ACAAACAGAGAAATATGACACACACACCGGTGTTGATAGTCAGGTGA'
        dinuc = seqsearch.MicrosatelliteMotif(motif.dinucleotide, 6, True, 'Dinucleotide')
        ac_repeat = dinuc.compiled[0]
        matches = ac_repeat.finditer(base_seq)
        i = 0
        for m in matches:
            i += 1
        assert i == 0
        
class Test_2_Sequence:
    def setUp(self):
        base_seq = Seq('TTGTGTAAAACTCGCCTGGAAGACTGATGAA')
        self.s = seqsearch.RegionSearch(base_seq)
    def test_2_forward(self):
        '''[RegionSearch] : ensure sequence is a proper FORWARD representation of itself'''
        assert self.s.forward == 'TTGTGTAAAACTCGCCTGGAAGACTGATGAA'
    def test_3_type(self):
        '''[RegionSearch] : ensure passed objects must be Seq instances'''
        base_seq = 'TTGTGTAAAACTCGCCTGGAAGACTGATGAA'
        tools.assert_raises(AssertionError, seqsearch.RegionSearch, base_seq)
    def tearDown(self):
        base_seq = None
        self.s = None

class Test_4_RegionSearch:
    def setUp(self):
        self.min_repeats = 6
        self.mononucleotide     = seqsearch.MicrosatelliteMotif(motif.mononucleotide, self.min_repeats, False, name='Mononucleotide')
        self.dinucleotide       = seqsearch.MicrosatelliteMotif(motif.dinucleotide, self.min_repeats, False, name='Dinucleotide')
        self.trinucleotide      = seqsearch.MicrosatelliteMotif(motif.trinucleotide, self.min_repeats, False, name = 'Trinucleotide')
        self.tetranucleotide    = seqsearch.MicrosatelliteMotif(motif.tetranucleotide, self.min_repeats, False, name = 'Tetranucleotide')
        self.pentanucleotide    = seqsearch.MicrosatelliteMotif(motif.pentanucleotide, self.min_repeats, False, name = 'Pentanucleotide')
        self.hexanucleotide     = seqsearch.MicrosatelliteMotif(motif.hexanucleotide, self.min_repeats, False, name = 'Hexanucleotide')
    
    def generic(self, msat_length, repeat_instance):
        handle = open('testRepeats.fa','rU')
        for seq in SeqIO.parse(handle, "fasta"):
            values = seq.id.split('_')
            motif_length, repeat_seq, expected_length, expected_span = int(values[0]), values[1], int(values[2]), eval(values[3])
            if motif_length == msat_length:
                search = seqsearch.RegionSearch(seq.seq)
                search.microsatellite(repeat_instance)
                if expected_length >= self.min_repeats:
                    observed_span = search.matches[repeat_seq][0]
                    observed_length = (observed_span[1]-observed_span[0])/msat_length
                    assert expected_span == observed_span, 'Mismatched repeat spans %s; expected: %s observed:%s' % (repeat_seq, expected_span, observed_span)
                    assert expected_length == observed_length, 'Mismatches repeat lengths %s, expected: %s, got: %s' % (repeat_seq, expected_length, observed_length)
                else:
                    assert search.matches == {}
                
    def test_1_simulated_perfect_repeats(self):
        '''[RegionSearch.microsatellite]:  searching for mononucs and checking length'''
        self.generic(1, self.mononucleotide)
        
    def test_2_simulated_perfect_repeats(self):
        '''[RegionSearch.microsatellite]:  searching for dinucs and checking length'''
        self.generic(2, self.dinucleotide)
    
    def test_3_simulated_perfect_repeats(self):
        '''[RegionSearch.microsatellite]:  searching for trinucs and checking length'''
        self.generic(3, self.trinucleotide)

    def test_4_simulated_perfect_repeats(self):
        '''[RegionSearch.microsatellite]:  searching for tetranucs and checking length'''
        self.generic(4, self.tetranucleotide)
    
    def test_5_simulated_perfect_repeats(self):
        '''[RegionSearch.microsatellite]:  searching for pentanucs and checking length'''
        self.generic(5, self.pentanucleotide)

    def test_6_simulated_perfect_repeats(self):
        '''[RegionSearch.microsatellite]:  searching for hexanucs and checking length'''
        self.generic(6, self.hexanucleotide)        
    