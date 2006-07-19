#!/usr/bin/env python

#-----------------------------------------------------------------------------------------------|
# MicrosatFinder - command line version for Python:  an update to N. Dean Pentcheff's           |
# Ephemeris 1.0 for Perl                                                                        |
# available from http://www.uga.edu/srel/DNA_Lab/ephemeris%201.0.bin.                           |
#                                                                                               |
# Copyright (C) 2005-2006 Brant C. Faircloth.  Modifications/recoding conducted solely by       |
# Brant C. Faircloth in March, 2005 porting program to Python (www.python.org) and recoding     |
# some parts.                                                                                   |                                                                                    
#                                                                                               |                                                                             
# This program is free software; you can redistribute it and/or modify it under the terms of the|
# GNU General Public License as published by the Free Software Foundation; either version 2 of  | 
# the License, or (at your option) any later version.                                           |
#                                                                                               |
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     |
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     |
# See the GNU General Public License for more details.                                          |
#                                                                                               |
# You should have received a copy of the GNU General Public License along with this program; if |
# not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     |
# 02111-1307, USA.                                                                              |
#-----------------------------------------------------------------------------------------------|

import os, string, sys, getopt, re

try:
    from Bio import Fasta
except:
    print "This program requires BioPython (http://biopython.org/) to be installed."
       
class progressBar:
    """ Creates a text-based progress bar. Call the object with the `print'
        command to see the progress bar, which looks something like this:
            
        [=======>        22%                  ]
        
        You may specify the progress bar's width, min and max values on init.
    """

    def __init__(self, minValue = 0, maxValue = 100, totalWidth=80):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, newAmount = 0):
        """ Update the progress bar with the new amount (with min and max
            values set at initialization; if it is over or under, it takes the
            min or max value as a default. """
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = "[>%s]" % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % ('='*allFull)
        else:
            self.progBar = "[%s>%s]" % ('='*(numHashes-1),
                                        ' '*(allFull-numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                                self.progBar[percentPlace+len(percentString):]
                                ])

    def __str__(self):
        return str(self.progBar)

    def __call__(self, value):
        """ Updates the amount, and writes to stdout. Prints a carriage return
            first, so it will overwrite the current line in stdout."""
        print '\r',
        self.updateAmount(value)
        sys.stdout.write(str(self))
        sys.stdout.flush()

class mods: 
    
    """Class representing DNA as a string sequence.""" 

    def complement(self, s):
            
            """return complementary dna sequence"""
            
            tab = string.maketrans('AGCTagct','TCGAtcga')
            output = string.translate(s, tab)
            return output
    
class search:
    
    def __init__(self):
         
        """Create DNA instance initialized to string s."""
        
        self.repeatUnits = {"mononucleotide":1, "dinucleotide":2, "trinucleotide":3, "tetranucleotide":4, "pentanucleotide":5, "hexanucleotide":6}
        self.minSize={"mononucleotide":'{9,}',"dinucleotide":'{6,}',"trinucleotide":'{4,}',"tetranucleotide":'{3,}', "pentanucleotide":'{3,}', "hexanucleotide":'{3,}'} #defines minimum size for repeat unit
        #defines repeat units (lowest alphabetical, unique, non-complementary) for which we are searching
        self.mononucleotide=['(A)','(C)'] 
        
        self.dinucleotide=['(AC)','(AG)',
                            '(AT)','(CG)'] 
        
        self.trinucleotide=['(AAC)','(AAG)','(AAT)',
                            '(ACC)','(ACG)','(ACT)',
                            '(AGC)','(AGG)','(ATC)',
                            '(CCG)']
                            
        self.tetranucleotide=['(AAAC)','(AAAG)','(AAAT)','(AACC)',
                            '(AACG)','(AACT)','(AAGC)','(AAGG)',
                            '(AAGT)','(ACAG)','(ACAT)','(ACCC)',
                            '(ACCG)','(ACCT)','(ACGC)','(ACGT)',
                            '(ACTC)','(ACTG)','(AGAT)','(AATC)',
                            '(AATG)','(AATT)','(ACGG)','(AGCC)',
                            '(AGCG)','(AGGC)','(AGGG)','(ATCC)',
                            '(ATCG)','(ATGC)','(CCCG)','(CCGG)']
        
        self.pentanucleotide=['(AAAAC)', '(AAAAG)', '(AAAAT)', '(AAACC)', '(AAACG)', 
                            '(AAACT)', '(AAAGC)', '(AAAGG)', '(AAAGT)', '(AAATC)', 
                            '(AAATG)', '(AAATT)', '(AACAC)', '(AACAG)', '(AACAT)', 
                            '(AACCC)', '(AACCG)', '(AACCT)', '(AACGC)', '(AACGG)',
                            '(AACGT)', '(AACTC)', '(AACTG)', '(AACTT)', '(AAGAC)',
                            '(AAGAG)', '(AAGAT)', '(AAGCC)', '(AAGCG)', '(AAGCT)',
                            '(AAGGC)', '(AAGGG)', '(AAGGT)', '(AAGTC)', '(AAGTG)', 
                            '(AATAC)', '(AATAG)', '(AATAT)', '(AATCC)', '(AATCG)', 
                            '(AATCT)', '(AATGC)', '(AATGG)', '(AATGT)', '(AATTC)', 
                            '(ACACC)', '(ACACG)', '(ACACT)', '(ACAGC)', '(ACAGG)', 
                            '(ACAGT)', '(ACATC)', '(ACATG)', '(ACCAG)', '(ACCAT)', 
                            '(ACCCC)', '(ACCCG)', '(ACCCT)', '(ACCGC)', '(ACCGG)', 
                            '(ACCGT)', '(ACCTC)', '(ACCTG)', '(ACGAG)', '(ACGAT)', 
                            '(ACGCC)', '(ACGCG)', '(ACGCT)', '(ACGGC)', '(ACGGG)', 
                            '(ACGTC)', '(ACTAG)', '(ACTAT)', '(ACTCC)', '(ACTCG)', 
                            '(ACTCT)', '(ACTGC)', '(ACTGG)', '(AGAGC)', '(AGAGG)', 
                            '(AGATC)', '(AGATG)', '(AGCAT)', '(AGCCC)', '(AGCCG)', 
                            '(AGCCT)', '(AGCGC)', '(AGCGG)', '(AGCTC)', '(AGGAT)', 
                            '(AGGCC)', '(AGGCG)', '(AGGGC)', '(AGGGG)', '(ATATC)', 
                            '(ATCCC)', '(ATCCG)', '(ATCGC)', '(ATGCC)', '(CCCCG)', 
                            '(CCCGG)', '(CCGCG)']
                            
        self.hexanucleotide=['(AAAAAC)', '(AAAAAG)', '(AAAAAT)', '(AAAACC)', '(AAAACG)', '(AAAACT)', 
                            '(AAAAGC)', '(AAAAGG)', '(AAAAGT)', '(AAAATC)', '(AAAATG)', '(AAAATT)',
                            '(AAACAC)', '(AAACAG)', '(AAACAT)', '(AAACCC)', '(AAACCG)', '(AAACCT)',
                            '(AAACGC)', '(AAACGG)', '(AAACGT)', '(AAACTC)', '(AAACTG)', '(AAACTT)', 
                            '(AAAGAC)', '(AAAGAG)', '(AAAGAT)', '(AAAGCC)', '(AAAGCG)', '(AAAGCT)', 
                            '(AAAGGC)', '(AAAGGG)', '(AAAGGT)', '(AAAGTC)', '(AAAGTG)', '(AAAGTT)', 
                            '(AAATAC)', '(AAATAG)', '(AAATAT)', '(AAATCC)', '(AAATCG)', '(AAATCT)', 
                            '(AAATGC)', '(AAATGG)', '(AAATGT)', '(AAATTC)', '(AAATTG)', '(AAATTT)', 
                            '(AACAAG)', '(AACAAT)', '(AACACC)', '(AACACG)', '(AACACT)', '(AACAGC)', 
                            '(AACAGG)', '(AACAGT)', '(AACATC)', '(AACATG)', '(AACATT)', '(AACCAC)', 
                            '(AACCAG)', '(AACCAT)', '(AACCCC)', '(AACCCG)', '(AACCCT)', '(AACCGC)', 
                            '(AACCGG)', '(AACCGT)', '(AACCTC)', '(AACCTG)', '(AACCTT)', '(AACGAC)', 
                            '(AACGAG)', '(AACGAT)', '(AACGCC)', '(AACGCG)', '(AACGCT)', '(AACGGC)', 
                            '(AACGGG)', '(AACGGT)', '(AACGTC)', '(AACGTG)', '(AACGTT)', '(AACTAC)', 
                            '(AACTAG)', '(AACTAT)', '(AACTCC)', '(AACTCG)', '(AACTCT)', '(AACTGC)', 
                            '(AACTGG)', '(AACTGT)', '(AACTTC)', '(AACTTG)', '(AAGAAT)', '(AAGACC)', 
                            '(AAGACG)', '(AAGACT)', '(AAGAGC)', '(AAGAGG)', '(AAGAGT)', '(AAGATC)', 
                            '(AAGATG)', '(AAGATT)', '(AAGCAC)', '(AAGCAG)', '(AAGCAT)', '(AAGCCC)', 
                            '(AAGCCG)', '(AAGCCT)', '(AAGCGC)', '(AAGCGG)', '(AAGCGT)', '(AAGCTC)', 
                            '(AAGCTG)', '(AAGCTT)', '(AAGGAC)', '(AAGGAG)', '(AAGGAT)', '(AAGGCC)', 
                            '(AAGGCG)', '(AAGGCT)', '(AAGGGC)', '(AAGGGG)', '(AAGGGT)', '(AAGGTC)', 
                            '(AAGGTG)', '(AAGTAC)', '(AAGTAG)', '(AAGTAT)', '(AAGTCC)', '(AAGTCG)', 
                            '(AAGTCT)', '(AAGTGC)', '(AAGTGG)', '(AAGTGT)', '(AATACC)', '(AATACG)', 
                            '(AATACT)', '(AATAGC)', '(AATAGG)', '(AATAGT)', '(AATATC)', '(AATATG)', 
                            '(AATATT)', '(AATCAC)', '(AATCAG)', '(AATCAT)', '(AATCCC)', '(AATCCG)', 
                            '(AATCCT)', '(AATCGC)', '(AATCGG)', '(AATCGT)', '(AATCTC)', '(AATCTG)', 
                            '(AATGAC)', '(AATGAG)', '(AATGAT)', '(AATGCC)', '(AATGCG)', '(AATGCT)', 
                            '(AATGGC)', '(AATGGG)', '(AATGGT)', '(AATGTC)', '(AATGTG)', '(AATTAC)', 
                            '(AATTAG)', '(AATTAT)', '(AATTCC)', '(AATTCG)', '(AATTGC)', '(ACACAG)', 
                            '(ACACAT)', '(ACACCC)', '(ACACCG)', '(ACACCT)', '(ACACGC)', '(ACACGG)', 
                            '(ACACGT)', '(ACACTC)', '(ACACTG)', '(ACAGAG)', '(ACAGAT)', '(ACAGCC)', 
                            '(ACAGCG)', '(ACAGCT)', '(ACAGGC)', '(ACAGGG)', '(ACAGGT)', '(ACAGTC)', 
                            '(ACAGTG)', '(ACATAG)', '(ACATAT)', '(ACATCC)', '(ACATCG)', '(ACATCT)', 
                            '(ACATGC)', '(ACATGG)', '(ACATGT)', '(ACCACG)', '(ACCACT)', '(ACCAGC)', 
                            '(ACCAGG)', '(ACCAGT)', '(ACCATC)', '(ACCATG)', '(ACCCAG)', '(ACCCAT)', 
                            '(ACCCCC)', '(ACCCCG)', '(ACCCCT)', '(ACCCGC)', '(ACCCGG)', '(ACCCGT)', 
                            '(ACCCTC)', '(ACCCTG)', '(ACCGAG)', '(ACCGAT)', '(ACCGCC)', '(ACCGCG)', 
                            '(ACCGCT)', '(ACCGGC)', '(ACCGGG)', '(ACCGGT)', '(ACCGTC)', '(ACCGTG)', 
                            '(ACCTAG)', '(ACCTAT)', '(ACCTCC)', '(ACCTCG)', '(ACCTCT)', '(ACCTGC)', 
                            '(ACCTGG)', '(ACGACT)', '(ACGAGC)', '(ACGAGG)', '(ACGAGT)', '(ACGATC)', 
                            '(ACGATG)', '(ACGCAG)', '(ACGCAT)', '(ACGCCC)', '(ACGCCG)', '(ACGCCT)', 
                            '(ACGCGC)', '(ACGCGG)', '(ACGCGT)', '(ACGCTC)', '(ACGCTG)', '(ACGGAG)', 
                            '(ACGGAT)', '(ACGGCC)', '(ACGGCG)', '(ACGGCT)', '(ACGGGC)', '(ACGGGG)', 
                            '(ACGTAG)', '(ACGTAT)', '(ACGTCC)', '(ACGTCG)', '(ACGTGC)', '(ACTAGC)', 
                            '(ACTAGG)', '(ACTAGT)', '(ACTATC)', '(ACTATG)', '(ACTCAG)', '(ACTCAT)', 
                            '(ACTCCC)', '(ACTCCG)', '(ACTCCT)', '(ACTCGC)', '(ACTCGG)', '(ACTCTC)', 
                            '(ACTCTG)', '(ACTGAG)', '(ACTGAT)', '(ACTGCC)', '(ACTGCG)', '(ACTGCT)', 
                            '(ACTGGC)', '(ACTGGG)', '(AGAGAT)', '(AGAGCC)', '(AGAGCG)', '(AGAGCT)', 
                            '(AGAGGC)', '(AGAGGG)', '(AGATAT)', '(AGATCC)', '(AGATCG)', '(AGATCT)', 
                            '(AGATGC)', '(AGATGG)', '(AGCAGG)', '(AGCATC)', '(AGCATG)', '(AGCCAT)', 
                            '(AGCCCC)', '(AGCCCG)', '(AGCCCT)', '(AGCCGC)', '(AGCCGG)', '(AGCCTC)', 
                            '(AGCCTG)', '(AGCGAT)', '(AGCGCC)', '(AGCGCG)', '(AGCGCT)', '(AGCGGC)', 
                            '(AGCGGG)', '(AGCTAT)', '(AGCTCC)', '(AGCTCG)', '(AGCTGC)', '(AGGATC)', 
                            '(AGGATG)', '(AGGCAT)', '(AGGCCC)', '(AGGCCG)', '(AGGCCT)', '(AGGCGC)', 
                            '(AGGCGG)', '(AGGGAT)', '(AGGGCC)', '(AGGGCG)', '(AGGGGC)', '(AGGGGG)', 
                            '(ATATCC)', '(ATATCG)', '(ATATGC)', '(ATCATG)', '(ATCCCC)', '(ATCCCG)', 
                            '(ATCCGC)', '(ATCCGG)', '(ATCGCC)', '(ATCGCG)', '(ATCGGC)', '(ATGCCC)', 
                            '(ATGCGC)', '(ATGGCC)', '(CCCCCG)', '(CCCCGG)', '(CCCGCG)', '(CCCGGG)', 
                            '(CCGCGG)', '(CCGGCG)']        

    def genericMethod(self, i, repeat):
        
        """generic method for finding various microsatellite repeats"""
        
        #print (("Finding repeats of size %s....\n") % (repeat))
        wildcard='[N]*' + i + '+'                                               # build wildcard string for N bases
        searchString = i + self.minSize[repeat] + wildcard                      # concatenates mononuc[i] and minimum size for repeat type
        compiledRegEx = re.compile(searchString, re.IGNORECASE)                 # compiles regex for concatenated values
        iterator = compiledRegEx.finditer(self.seq)                             # sets up iterative repeat finder
        compIterator = compiledRegEx.finditer(mods().complement(self.seq))
        for match in iterator:
            bases = match.span()                                                # give start/end bases of repeat
            length = (bases[1] - bases[0]) / self.repeatUnits[repeat]           # determine number of repeats for given msat type
            self.msatResults[bases[0]+1] = ('%s repeat %s^%s found between bases %s and %s.') % (repeat, i, length, bases[0]+1, bases[1]+1)
        for match in compIterator:                                              # do the same on the reverse complement of the sequence
            bases = match.span()
            length = (bases[1] - bases[0]) / self.repeatUnits[repeat]
            seq = match.group()
            self.msatResults[bases[0]+1] = ('Reverse complement of %s repeat %s, %s^%s found between bases %s and %s.') % (repeat, i, mods().complement(i), length, bases[0]+1, bases[1]+1)

    def ephemeris(self, s, type):
        
        """Searches for microsatellite sequences (mononucleotide, dinucleotide, trinucleotide, tetranucleotide) in DNA string"""        
        
        self.seq = s
        self.msatResults={}                                 # we will store output for each repeat in dictionary keyed on start base #
        
        if type == 'tetra':
            for repeatClass in ['self.mononucleotide','self.dinucleotide', 'self.trinucleotide', 'self.tetranucleotide']:
                for i in repeatClass:
                    self.genericMethod(i,repeatClass.lstrip('self.'))        
            #for i in self.dinucleotide:
            #    self.genericMethod(i,"dinucleotide")
            #for i in self.trinucleotide:
            #    self.genericMethod(i,"trinucleotide")
            #for i in self.tetranucleotide:
            #    self.genericMethod(i, "tetranucleotide")
        elif type == 'penta':
            for repeatClass in ['self.mononucleotide','self.dinucleotide', 'self.trinucleotide', 'self.tetranucleotide', 'self.pentanucleotide']:
                for i in repeatClass:
                    self.genericMethod(i,repeatClass.lstrip('self.'))
            #for i in self.mononucleotide:
            #    self.genericMethod(i,"mononucleotide")        
            #for i in self.dinucleotide:
            #    self.genericMethod(i,"dinucleotide")
            #for i in self.trinucleotide:
            #    self.genericMethod(i,"trinucleotide")
            #for i in self.tetranucleotide:
            #   self.genericMethod(i, "tetranucleotide")
            #for i in self.pentanucleotide:
            #    self.genericMethod(i, "pentanucleotide")
        else:
            for repeatClass in ['self.mononucleotide','self.dinucleotide', 'self.trinucleotide', 'self.tetranucleotide', 'self.pentanucleotide', self.hexanucleotide]:
                for i in repeatClass:
                    self.genericMethod(i,repeatClass.lstrip('self.'))
            #for i in self.mononucleotide:
            #    self.genericMethod(i,"mononucleotide")        
            #for i in self.dinucleotide:
            #    self.genericMethod(i,"dinucleotide")
            #for i in self.trinucleotide:
            #    self.genericMethod(i,"trinucleotide")
            #for i in self.tetranucleotide:
            #   self.genericMethod(i, "tetranucleotide")
            #for i in self.pentanucleotide:
            #    self.genericMethod(i, "pentanucleotide")
            #for i in self.hexanucleotide:
            #    self.genericMethod(i, "hexanucleotide")

        return self.msatResults

def Usage():
    print "microsatFinder [-i] 'input filename' [-o] 'output filename' [-s] 'search type (tetra | penta | hexa)' [-h] help"
    sys.exit()

def getUserFiles():
    optlist, list = getopt.getopt(sys.argv[1:], 'i:o:s:vh')
    output = ''                                                     # set output to empty
    searchType = 'hexa'                                             #set default search to hexanucleotide (e.g. all)
    i=0
    if optlist:
        for opt in optlist:
            if optlist[i][0] == '-i':
                inFile=optlist[i][1]
            elif optlist[i][0] == '-o':
                outFile=optlist[i][1]
            elif optlist[i][0] == '-s':
                searchType=optlist[i][1]
            elif optlist[i][0] == '-h':
                Usage()
            i+=1
        try:
            inFile = os.path.abspath(string.strip(inFile))        # have to strip whitespace characters for dragging folders
            #fileList=getFiles(input)                            # calls function above
        except:
            print 'File/Directory does not exist!'
            sys.exit()
        if searchType in ['tetra','penta','hexa']:
                print (('\nYou are searching for all %sNUCLEOTIDE and smaller repeats.') % (searchType.upper()))
        else:
            print "Please choose 'tetra'|'penta'|'hexa' or leave blank for the default (hexa).\n"
            sys.exit()
        if not output:
            outFile=os.path.join(os.path.dirname(inFile),'output.txt')
        else:
            try:
                os.path.isdir(os.path.dirname(os.path.abspath(outFile)))
            except:
                print 'This is not a valid path.  Make sure you have entered the path correctly.'
                sys.exit()
        print ('\nOutput file written to %s.\n') % (os.path.abspath(outFile))
    
    else:
        print "Please enter options on the command line in the form of\n\n./msatCommandFasta.py -i 'infile.txt' -o 'outfile.txt'\n\nTry ./msatCommand -h for help."
    
    
    return inFile, outFile, searchType
    
def readInfo(inFile, outFile, repeatChoice):
    
    print 'Scanning sequences for microsatellites.  Percent Complete:'
    
    file = open(inFile)
    data = file.readlines()
    file.close()
    length = len(data)/2
    interval = 1./length * 100.
    
    prog = progressBar(0,100,80)
    
    file=open(outFile,'w')                                   # opens file for output - append only to keep from overwriting
    file.write('Microsatellite repeats found in the following sequences: \n\n')

    parser = Fasta.RecordParser()
    infile = open(inFile)
    iterator = Fasta.Iterator(infile, parser)
    i = 0
    while 1:
        record = iterator.next()
        if not record:
            break
            infile.close()
            file.close()
        dataOut=search().ephemeris(record.sequence, repeatChoice) 
        dictKeys=dataOut.keys()
        dictKeys.sort()                                     # sorts keys so bp locations will be in order
        if dictKeys:
            #file.write(('%s>>%s %s') % ('\n', record.title, '\n'))
            for k in dictKeys:                                  # writes dict values for sorted keys to output file
                dataList = dataOut[k].split()
                file.write(('%s\t%s\t%s\t%s\n') % (record.title, ' '.join(dataList[:-7]), dataList[-7], ' '.join(dataList[-6:])))
            file.write(('---------------------------------------%s') % ('\n'))
        i += interval
        prog(i)

if __name__ == '__main__':
    userInput,userOutput,userRepeatChoice = getUserFiles()
    readInfo(userInput, userOutput, userRepeatChoice)
    print '\n'