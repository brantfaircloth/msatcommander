#!/usr/bin/env python
# encoding: utf-8
"""
msatCommander.py

Copyright (c) 2005-2007 Brant C. Faircloth. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation; either version 2 of   
the License, or (at your option) any later version.                                           
                                                                                           
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     
See the GNU General Public License for more details.                                          
                                                                                           
You should have received a copy of the GNU General Public License along with this program; if 
not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     
02111-1307, USA.
 
"""

import string, re, time, csv, os.path, getopt, sys
from Bio.SeqIO import SequenceIterator
       
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

class search:
    def __init__(self,s):
        self.seq = s
    
    def revComplement(self, s):
        """return complementary dna sequence"""
        bases = string.maketrans('AGCTagct','TCGAtcga')
        # translate it
        output = string.translate(s, bases)
        # reverse it
        output = output[::-1]
        return output
    
    def genericMethod(self, passedSearchClass, repeat):
        """generic method for finding various microsatellite repeats"""
        # quick and dirty means of indexing position in passedSearchClass
        compiledRegExPos = 0                                                   
        sLength = len(self.seq)
        for compiledRegEx in passedSearchClass:
            iterator = compiledRegEx.finditer(self.seq)
            compIterator = compiledRegEx.finditer(self.revComplement(self.seq))
            name = repeat + 'Letters'
            i = allRepeatClasses[name][compiledRegExPos]
            # list to help remove duplicates
            baseList = []
            # look for matches in forward sequence
            for match in iterator:
                bases = match.span()                                                
                # give start/end bases of repeat
                length = (bases[1] - bases[0]) / repeatUnits[repeat]
                self.msatResults[bases[0]+1] = ('%s repeat %s^%s found between bases %s and %s.') % (repeat, i, length, bases[0]+1, bases[1]+1)
                baseList.append(bases[0])
            # do the same on the reverse complement of the sequence
            for match in compIterator:                                              
                bases = match.span()
                # add if statement to remove repetitive finds...
                if sLength - bases[1] not in baseList:
                    length = (bases[1] - bases[0]) / repeatUnits[repeat]
                    #seq = match.group()            
                    i = i.strip('()')
                    self.msatResults[bases[0]+1] = ('Reverse complement of %s repeat %s, (%s)^%s found between bases %s and %s.') % (repeat, i, self.revComplement(i), length, sLength - bases[1]+1, sLength - bases[0]+1)
            compiledRegExPos+=1

    def ephemeris(self, type):
        
        """Searches for microsatellite sequences (mononucleotide, dinucleotide, 
        trinucleotide, tetranucleotide) in DNA string"""
        
        # we will store output for each repeat in dictionary keyed on start base #
        self.msatResults={}                                 
        #print 's.mono ', self.mononucleotide
        # moved All (slow!) here to keep from double searching
        if 'All' in type:                          
            self.genericMethod(mononucleotide, 'mononucleotide')
            self.genericMethod(dinucleotide,'dinucleotide')
            self.genericMethod(trinucleotide,'trinucleotide')
            self.genericMethod(tetranucleotide,'tetranucleotide')
            self.genericMethod(pentanucleotide,'pentanucleotide')
            self.genericMethod(hexanucleotide,'hexanucleotide')      
        else:
            for searchClass in type:
                if str(searchClass) == 'Mononucleotide':
                    self.genericMethod(mononucleotide, 'mononucleotide')        
                elif str(searchClass) == 'Dinucleotide':
                    self.genericMethod(dinucleotide,'dinucleotide')              
                elif str(searchClass) == 'Trinucleotide':
                    self.genericMethod(trinucleotide,'trinucleotide')
                elif str(searchClass) == 'Tetranucleotide':             
                    self.genericMethod(tetranucleotide,'tetranucleotide')
                elif str(searchClass) == 'Pentanucleotide':               
                    self.genericMethod(pentanucleotide,'pentanucleotide')
                elif str(searchClass) == 'Hexanucleotide':
                    self.genericMethod(hexanucleotide,'hexanucleotide')
        return self.msatResults

class excelSingleSpace:
    """class for csv module to work correctly"""
    delimiter = ','
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r'

class fileFunctions:
    
    def __init__(self):
        # initialize to defaults
        self.outfile, self.selection = '','All'
        # get options passed on CL
        try:
            opts, args = getopt.getopt(sys.argv[1:], 'i:o:s:vh')
        except getopt.GetoptError:
            self.Usage()
            sys.exit(2)
        for o,a in opts:
            if o in ("-i","--input"):self.infile = a
            if o in ("-o","--output"):self.outfile = a
            if o in ("-s","--search"):self.selection = a
            if o in ("-h","--help"):
                self.Usage()
                sys.exit()
    
    def Usage(self):
        print """
        Usage:
        ==========
        microsatFinder [-i|--input] FILE [-o|--output] FILE [-s|--search] TYPE [-h|--help]
        """
        sys.exit()
    
    def fileExceptions(self):
        print "\n***********************************"
        print "msatcommander\ncommand-lineversion 0.4.5"
        print "**********************************"
        try:
            # have to strip whitespace characters for dragging folders
            self.infile = os.path.abspath(string.strip(self.infile))
            # calls function above
            #fileList=getFiles(input)
        except:
            print 'File/Directory does not exist!\n\nExiting.'
            sys.exit()
        if self.selection in ['tetra','penta','hexa', 'All']:
                print (('\nYou are searching for %s microsatellite repeats.\n') % (self.selection.upper()))
        else:
            print "Please choose 'tetra'|'penta'|'hexa' or leave blank for the default (hexa).\n"
            sys.exit()
        if not self.outfile:
            self.outfile=os.path.join(os.path.dirname(self.infile),'output.txt')
            print (("Your file will be saved as %s\n") % (self.outfile))
            userCheck = raw_input("Are you sure (Y/n)? ")
            if userCheck not in ["y","Y"]:
                print "\nDoing nothing. Exiting."
                sys.exit()
        else:
            try:
                os.path.isdir(os.path.dirname(os.path.abspath(self.outfile)))
            except:
                print 'This is not a valid path.  Exiting.'
                sys.exit()
    
    def printRunData(self):
        for choice in self.selection:
            userChoice = (('%s ') % (choice))
        userChoice = (('You searched for all %s microsatellite repeats in the above sequences') % (userChoice))
        runTime = time.time() - self.startTime
        runTime = (('Time for execution = %f sec') % (runTime))
        self.csvWriter.writerows([
            [''],
            [userChoice],
            ["Sequences containing repeats:", sum(self.overallRepeatSequences)],
            ["TOTAL number of repeats found:", sum(self.overallRepeats)],
            ["Sequences searched for repeats:", self.sequenceCount],
            [''],
            [runTime]
            ])

    def RunSearch(self):
        print 'Scanning sequences for microsatellites.  Percent Complete:'

        handle = open(self.infile,"rU")
        length = 0
        for record in SequenceIterator(handle,"fasta"):
            length += 1 
        handle.close()
        #length = len(data)/2
        interval = 1./length * 100.

        prog = progressBar(0,100,80)
        
        self.startTime = time.time()
        # opens file for output (overwrite existing)
        file=open(self.outfile,'w')                                   
        self.csvWriter = csv.writer(file, dialect = excelSingleSpace)
        self.csvWriter.writerow(['Clone','Repeat Info','Repeat Count','Location','Start BP','','End BP'])
        # new-style biopython SeqIO iterator
        handle = open(self.infile,"rU")
        # initialize variable for number of sequence searched
        self.sequenceCount, startInt = 0,0
        # initialize list for sequence containing repeats
        self.overallRepeatSequences = []
        # overall repeats list
        self.overallRepeats = []
        for record in SequenceIterator(handle,"fasta"):
            sequence = search(record.seq.tostring())
            dataOut=sequence.ephemeris(self.selection) 
            dictKeys=dataOut.keys()
            # sorts keys so bp locations will be in order by clone name
            dictKeys.sort()                                    
            if dictKeys:
                # add item for sequence containing repeat
                self.overallRepeatSequences.append(1)
                # writes dict values for sorted keys to output file
                for k in dictKeys:                                  
                    self.overallRepeats.append(1)
                    dataList = dataOut[k].split()
                    cloneResults = [record.id, ' '.join(dataList[:-7]), dataList[-7], ' '.join(dataList[-6:-3]), dataList[-3], dataList[-2], dataList[-1]]
                    self.csvWriter.writerow(cloneResults)
            elif not dictKeys and self.noRepeats:
                cloneResults = [record.id, "No repeats found"]        
                self.csvWriter.writerow(cloneResults)
            self.sequenceCount +=1
            startInt += interval
            prog(startInt)       
        self.printRunData()
        print ('\n\nOutput file written to:\n%s') % (os.path.abspath(self.outfile))
        # close open files
        handle.close()
        file.close()

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
    #------------------------------------------------------------------------
    # number of units (letters) in each repeat class - used to determine size
    #------------------------------------------------------------------------
    repeatUnits = {"mononucleotide":1, "dinucleotide":2, "trinucleotide":3, "tetranucleotide":4, "pentanucleotide":5, "hexanucleotide":6}
    #------------------------------------------------------------------------
    # minimum number of repeats to qualify an area as being a microsatellite.
    # Number is one less than actual value
    #------------------------------------------------------------------------
    minSize={"mononucleotide":'{9,}',
    "dinucleotide":'{6,}',
    "trinucleotide":'{4,}',
    "tetranucleotide":'{3,}', 
    "pentanucleotide":'{3,}', 
    "hexanucleotide":'{3,}'}
    #---------------------------------------------------------------------------------------------
    # All repeat units (lowest alphabetical, unique, non-complementary) for which we are searching 
    # these are compiled into regular expressions below...
    #---------------------------------------------------------------------------------------------
    allRepeatClasses = {'mononucleotideLetters':['(A)','(C)'], 
    'dinucleotideLetters':['(AC)','(AG)',
                        '(AT)','(CG)'], 
    'trinucleotideLetters':['(AAC)','(AAG)','(AAT)',
                        '(ACC)','(ACG)','(ACT)',
                        '(AGC)','(AGG)','(ATC)',
                        '(CCG)'],
    'tetranucleotideLetters':['(AAAC)','(AAAG)','(AAAT)','(AACC)',
                        '(AACG)','(AACT)','(AAGC)','(AAGG)',
                        '(AAGT)','(ACAG)','(ACAT)','(ACCC)',
                        '(ACCG)','(ACCT)','(ACGC)','(ACGT)',
                        '(ACTC)','(ACTG)','(AGAT)','(AATC)',
                        '(AATG)','(AATT)','(ACGG)','(AGCC)',
                        '(AGCG)','(AGGC)','(AGGG)','(ATCC)',
                        '(ATCG)','(ATGC)','(CCCG)','(CCGG)'],
    'pentanucleotideLetters':['(AAAAC)', '(AAAAG)', '(AAAAT)', '(AAACC)', '(AAACG)', 
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
                        '(CCCGG)', '(CCGCG)'],
    'hexanucleotideLetters':['(AAAAAC)', '(AAAAAG)', '(AAAAAT)', '(AAAACC)', '(AAAACG)', '(AAAACT)', 
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
                        '(CCGCGG)', '(CCGGCG)']} # closing squiggly for Letters dictionary
    #-----------------------------------------------------------------------------------------------------------------------------------
    # here:  go ahead and compile all regular expressions ONCE instead of at each pass (this is major part of speedup 2X mentiond above)
    #-----------------------------------------------------------------------------------------------------------------------------------
    for repeatClass in allRepeatClasses.keys():
        tempList = []
        for repeatType in allRepeatClasses[repeatClass]: 
            wildcard='[N]*' + repeatType + '+'                                               # build wildcard string for N bases
            searchString = repeatType + minSize[repeatClass[:-7]] + wildcard            # concatenates mononuc[i] and minimum size for repeat type
            compiledRegEx = re.compile(searchString, re.IGNORECASE)
            tempList.append(compiledRegEx)
        if repeatClass[:-7]=='mononucleotide':
            mononucleotide = tempList
        elif repeatClass[:-7]=='dinucleotide':
            dinucleotide = tempList
        elif repeatClass[:-7]=='trinucleotide':
            trinucleotide = tempList
        elif repeatClass[:-7]=='tetranucleotide':
            tetranucleotide = tempList
        elif repeatClass[:-7]=='pentanucleotide':
            pentanucleotide = tempList
        elif repeatClass[:-7]=='hexanucleotide':
            hexanucleotide = tempList
    #-----------------------------------------------------------------------
    # end global variable definition and compilation and begin main GUI loop    
    #-----------------------------------------------------------------------
    userData = fileFunctions()
    userData.fileExceptions()
    userData.RunSearch()
    #main()
    #userInput,userOutput,userRepeatChoice = getUserFiles()
    #readInfo(userInput, userOutput, userRepeatChoice)
    #print '\n'