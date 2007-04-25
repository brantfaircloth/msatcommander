#!/usr/bin/env python
# encoding: utf-8
"""
mscprimertag.py

msatCommander :: Copyright (c) 2005-2007 Brant C. Faircloth. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation; either version 2 of   
the License, or (at your option) any later version.                                           
                                                                                           
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     
See the GNU General Public License for more details.                                          
                                                                                           
You should have received a copy of the GNU General Public License along with this program; if 
not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA     
02111-1307, USA.

This program makes use of primer3, which is released under the *NEW* BSD License, 
here are the terms:

==================

Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky
All rights reserved.

Most of primer3 is released under the following _new_ BSD license:

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
   * Neither the names of the copyright holders nor contributors may
be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The oligtm library and tests are released under the GPL.  See
file src/gpl.txt or go to http://www.gnu.org/licenses/gpl.txt.
 
==================

"""
import sys, os, string, csv, subprocess, time

class mods: 
    """class defining several methods we will use to modify dna sequence.""" 
 
    def __init__(self, s): 
        """create DNA instance initialized to string s.""" 
        self.seq = s 
        
    def complement(self):
        """returns complementary dna sequence"""
        #--------------------------------------------------------------------------------------
        # first, we capitalize the sequence (if lowercase), then we translate it to its reverse
        # complement.
        #--------------------------------------------------------------------------------------
        if self.seq.islower():self.seq = self.seq.upper() 
        tab = string.maketrans('AGCTN','TCGAN')
        output = string.translate(self.seq, tab)
        return output

    def reverseComplement(self):
        """returns reverse complementary dna sequence"""
        #-----------------------------------------------------------------------------------------
        # the self.seq (sequence instance) must be converted to a list, reversed, put back together
        # as a string, and then complemented.
        #-----------------------------------------------------------------------------------------
        listS = list(self.seq)
        listS.reverse()
        self.seq = (''.join(listS))
        return self.complement()
        
class tagMain:
    """main class for sequence building and input file construction"""
    #--------------------------------------------------------------------------------------
    # This class is the main class we will use for building the temporary outfile for use
    # with primer3.  It is subclassed by several other classes and that allows us to only
    # define the main function 1X, thus saving us some code.
    #--------------------------------------------------------------------------------------
    def __init__(self, tag, upper, lower):
        self.fragTag=tag
        #---------------------------------------------------------------------------------
        # The vector sequence used here is the first 250 bp of Enterobacteria phage lambda
        # GenBank gi:9626243.  It really serves no other purpose than a spacer to put
        # between the primers we will be testing.
        #---------------------------------------------------------------------------------
        self.lamVector='GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGGCTTTTTGGCCTCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCAGCTGGCTGACATTTTCGGTGCGAGTATCCGTACCATTCAG'
        self.upperP = upper
        self.lowerP = lower
        self.lowerPRevComp = mods(self.lowerP).reverseComplement()
        #-------------------------------------------------------------------------------------
        # The above is essentially setting the variables which will be used in the subclasses.
        # That way, we don't have to define them again for each new function
        #-------------------------------------------------------------------------------------
        
    def commonBases(self, tag, primer):
        """returns tag; checks primer and tag for common bases and removes them"""
        self.tag = tag
        self.common = False
        i=5
        while i > 0:
            if self.tag.endswith(primer[0:i]):
                self.tag = self.tag.rstrip(primer[0:i])
                self.common = True
                break
            else:
                i -= 1
        return self.tag 
    
    def writeFile(self, outDir):
        """function builds temporary output file to be used by primer 3, subclassed below"""
        self.tempPrimer3TagFile = outDir + 'primer3_tag_temp.txt'
        self.outFile=open(self.tempPrimer3TagFile,'w')
        self.outFile.write(('SEQUENCE=%s\n') % (self.seqBuild()))
        if self.common:
            self.outFile.write(('PRIMER_COMMENT=Common Bases Modified, Tag %s\n') % (self.tag))
        self.tagAdd()
        #---------------------------------------------------------------------------------
        # Set (or unset) many default values/weights for primer3.  Range must be set high, 
        # TMs are set ridiculously high as the tag adds a lot of temp.
        #---------------------------------------------------------------------------------
        features = ['PRIMER_PRODUCT_SIZE_RANGE=100-500',
        'PRIMER_MIN_TM=45.0',
        'PRIMER_MAX_TM=100.0',
        'PRIMER_SELF_ANY=6.0',
        'PRIMER_SELF_END=2.5',
        'PRIMER_MIN_SIZE=15',
        'PRIMER_MAX_SIZE=45',
        'PRIMER_PICK_ANYWAY=1',
        'PRIMER_EXPLAIN_FLAG=1',
        'PRIMER_WT_TM_GT=0.0',
        'PRIMER_WT_TM_LT=0.0',
        'PRIMER_WT_SIZE_GT=0.0',
        'PRIMER_WT_SIZE_LT=0.0',
        'PRIMER_WT_COMPL_ANY=1.0',
        'PRIMER_WT_COMPL_END=1.0',
        'PRIMER_WT_END_STABILITY=1.0',
        'PRIMER_PAIR_WT_COMPL_ANY=1.0',
        'PRIMER_PAIR_WT_COMPL_END=1.0',
        'PRIMER_WT_GC_PERCENT_LT=0.0',   #do we want to keep GC?
        'PRIMER_WT_GC_PERCENT_LT=0.0',   #do we want to keep GC? (do we want to define GC?)
        '=']
        for feat in features:
            self.outFile.write(('%s\n') % (feat))
        self.outFile.close()
        #------------------------
        # End primer3 file values
        #------------------------

class tagFront(tagMain):
    """subclass of tagMain; front tag (forward/upper primer) sequence building and input file creation"""
    
    def seqBuild(self):
        """builds contig for primer selection"""
        self.fragTag = self.commonBases(self.fragTag, self.upperP)
        return self.fragTag + self.upperP + self.lamVector + self.lowerPRevComp
        #------------------------------------------------------------------------------------------------------
        # we are 'building' a contig for primer3 that contains the tag we want to use, the upper primer we are 
        # testing, the vector dna, and the reverse complement of our lower primer.
        #------------------------------------------------------------------------------------------------------
    
    def tagAdd(self):
        """builds upper primer with tag and sets lower primer value"""
        self.outFile.write(('PRIMER_LEFT_INPUT=%s\n') % (self.fragTag + self.upperP))
        self.outFile.write(('PRIMER_RIGHT_INPUT=%s\n') % (self.lowerP))
        return 

class tagRear(tagMain):
    """subclass of tagMain; rear tag (reverse/lower primer) sequence building and input file creation"""
    
    def seqBuild(self):
        """builds contig for primer selection"""
        self.fragTag = self.commonBases(self.fragTag, self.lowerP)
        self.fragTagRevComp = mods(self.fragTag).reverseComplement()
        #----------------------------------------------------------------------------
        # the above is required so we can generate the reverseComp of fragTag *after*
        # we remove common bases (before would hose us)
        #----------------------------------------------------------------------------
        return self.upperP + self.lamVector + self.lowerPRevComp + self.fragTagRevComp
    
    def tagAdd(self):
        """builds lower primer with tag and sets upper primer value"""
        self.outFile.write(('PRIMER_LEFT_INPUT=%s\n') % (self.upperP))
        self.outFile.write(('PRIMER_RIGHT_INPUT=%s\n') % (self.fragTag + self.lowerP))
        return

class primerInfo:
    """this class will test out tagged primers"""    
    def primer3Run(self, tempfile):
        if os.name == 'posix':
            pathToPrimer3 = (('./primer3_core < \"%s\"') % (tempfile))
        elif os.name == 'nt':
            pathToPrimer3 = (('primer3_core < \"%s\"') % (tempfile))
        p=subprocess.Popen(pathToPrimer3,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True).stdout
        primer3Output = p.read()
        p.close()
        primer3Output = primer3Output.split('\n')
        #--------------------------------------------------------------------------------------
        # here we call primer3, read the stdout from primer3, and put that output into a
        # string.  The chop up the string and put it in a dictionary - ahhhhh, dictionaries.
        #--------------------------------------------------------------------------------------
        primerDataDict={}
        for item in primer3Output[0:len(primer3Output)-2]:      #drop last 2 lines - this is an '=' and a space
            splitItem=string.split(item,'=')
            try:
                splitItem[1]=float(splitItem[1])                #convert what we can to floats
            except:
                pass
            primerDataDict[splitItem[0]]=splitItem[1]     
        stuffToKeep=['PRIMER_LEFT_END_STABILITY', 'PRIMER_LEFT_INPUT', 'PRIMER_LEFT_PENALTY', 'PRIMER_LEFT_SELF_ANY','PRIMER_LEFT_SELF_END', 'PRIMER_PAIR_COMPL_ANY', 'PRIMER_PAIR_COMPL_END', 'PRIMER_PAIR_PENALTY', 'PRIMER_RIGHT_END_STABILITY', 'PRIMER_RIGHT_INPUT', 'PRIMER_RIGHT_PENALTY', 'PRIMER_RIGHT_SELF_ANY', 'PRIMER_RIGHT_SELF_END', 'PRIMER_WARNING','PRIMER_COMMENT']
        for item in primerDataDict.keys():
            if item not in stuffToKeep:
                try:
                    del primerDataDict[item]
                except:
                    pass
        try:
            os.remove(tempfile)
        except OSError:
            time.sleep(0.1)
            os.remove(tempfile)
        return primerDataDict

class excelSingleSpace:
    """class for csv module to work correctly"""
    delimiter = ','
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r'

class primer3Tag:
    def __init__(self, outDir, infile):
        self.outDir = outDir
        self.primerFile = infile
        if os.name == 'posix':
            self.tagResults = self.outDir + '/primer3/primerTag.csv'
        elif os.name == 'nt':
            self.tagResults = self.outDir + '\\primer3\\primerTag.csv'
    def getData(self):
        """get the freaking data"""
        dataDialect = csv.Sniffer().sniff(self.primerFile)          # try to determine dialect automagically
        #headerTrue = csv.Sniffer().has_header(dataFile)            # try to determine if file has header row
        primerData = open(self.primerFile,'rU')                     # open datafile for read
        primerReader = csv.reader(primerData, delimiter=',', dialect = dataDialect)  # setup csv.reader
        # we know there is a header, skip it...
        primerReader.next()
        # keep to list here as clone names may not be unique
        self.primerData = []
        for row in primerReader:
            clone,left,right = row[0],row[1],row[5]
            self.primerData.append((clone,left,right))
        primerData.close()

    def writeData(self, clone, primer, header, output):
        if 'PRIMER_PAIR_PENALTY' in primer.keys():
            #pdb.set_trace()
            buildList = [clone]
            for element in header[1:]:
                try:
                    try:
                        elementFloat = float(primer[element])
                    except:
                        elementFloat = primer[element]
                    buildList.append(elementFloat)
                except:
                    buildList.append('')
            #if common:
            #    buildList.append('Common bases modified')
            else: buildList.append('')
            #print buildList
            output.writerow(buildList)
        else: pass

    def run(self, cag, m13, custom):
        outputFile = self.tagResults
        outputData = open(outputFile,'w')                              # open datafile for read
        output_writer = csv.writer(outputData, dialect = excelSingleSpace)
        header = [
        'Clone',
        'PRIMER_LEFT_INPUT',
        'PRIMER_LEFT_SELF_ANY',
        'PRIMER_LEFT_SELF_END',
        'PRIMER_RIGHT_INPUT',
        'PRIMER_RIGHT_SELF_ANY',
        'PRIMER_RIGHT_SELF_END',
        'PRIMER_PAIR_COMPL_ANY',
        'PRIMER_PAIR_COMPL_END',
        'PRIMER_WARNING',
        'PRIMER_COMMENT']
        output_writer.writerow(header)

        for item in self.primerData:
            if item[0] == '' or item[1]=='':
                pass
            else:
                #print item[0]
                resultsDict={}
                #print item[0], item[1], item[2]
                if cag:
                    cagTag = 'CAGTCGGGCGTCATCA'
                    tCagF = tagFront(cagTag, item[1], item[2])
                    tCagF.writeFile(self.outDir)
                    resultsDict[0]=primerInfo().primer3Run(tCagF.tempPrimer3TagFile)
                    tCagR = tagRear(cagTag, item[1], item[2])
                    tCagR.writeFile(self.outDir)
                    resultsDict[1]=primerInfo().primer3Run(tCagR.tempPrimer3TagFile)
                if m13:
                    m13rTag = 'GGAAACAGCTATGACCAT'
                    tM13F = tagFront(m13rTag, item[1], item[2])
                    tM13F.writeFile(self.outDir)
                    resultsDict[2]=primerInfo().primer3Run(tM13F.tempPrimer3TagFile)
                    tM13R = tagRear(m13rTag, item[1], item[2])
                    tM13R.writeFile(self.outDir)
                    resultsDict[3]=primerInfo().primer3Run(tM13R.tempPrimer3TagFile)
                if custom:
                    custom = custom.upper()
                    tCustom = tagFront(custom, item[1], item[2])
                    tCustom.writeFile(self.outDir)
                    resultsDict[4]=primerInfo().primer3Run(tCustom.tempPrimer3TagFile)
                    tCustomR = tagRear(custom, item[1], item[2])
                    tCustomR.writeFile(self.outDir)
                    resultsDict[5]=primerInfo().primer3Run(tCustomR.tempPrimer3TagFile)
                primerPicks=[]
                for element in resultsDict:
                    try:
                        primerPicks.append((resultsDict[element]['PRIMER_PAIR_PENALTY'],element))                                             
                    except:
                        primerPicks.append(('1000',1000))  #insert arbitrarily bad value to list for primer w/ no results
                primerPicks.sort()
                if primerPicks[0][1] != 1000:
                    goodPrimer = resultsDict[primerPicks[0][1]]        #choose MINIMUM primer penalty value and keep that record     
                    self.writeData(item[0],goodPrimer,header,output_writer)
                else:
                    pass
        outputData.close()                                     
