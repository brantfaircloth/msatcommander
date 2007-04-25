#!/usr/bin/env python
# encoding: utf-8
"""
mscfunc.py

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
 
"""

import string, re, wx, time, csv, repeatClasses
#from Bio import Fasta
from Bio import SeqIO

class repeats:
    def createSearchPattern(self, mLen, dLen, tLen, ttLen, pLen, hLen):
        # number of units (letters) in each repeat class - used to determine size
        self.units = {"mononucleotide":1,
        "dinucleotide":2,
        "trinucleotide":3,
        "tetranucleotide":4,
        "pentanucleotide":5,
        "hexanucleotide":6}
        # setup length in Dict for regex
        mLen = (('{%s,}') % (int(mLen) - 1))
        dLen = (('{%s,}') % (int(dLen) - 1))
        tLen = (('{%s,}') % (int(tLen) - 1))
        ttLen = (('{%s,}') % (int(ttLen) - 1))
        pLen = (('{%s,}') % (int(pLen) - 1))
        hLen = (('{%s,}') % (int(hLen) - 1))
        # minimum number of repeats to qualify an area as being a microsatellite.
        # number is one less than actual value
        minSize={"mononucleotide":mLen,
        "dinucleotide":dLen,
        "trinucleotide":tLen,
        "tetranucleotide":ttLen, 
        "pentanucleotide":pLen, 
        "hexanucleotide":hLen}
        # All repeat units (lowest alphabetical, unique, non-complementary) for which we are searching 
        # these are compiled into regular expressions below...
        #--------------------------------------------------
        # here:  go ahead and compile all regular expressions ONCE instead of at each pass 
        # (this is major part of speedup 2X mentiond above)
        for repeatClass in repeatClasses.allRepeats.keys():
            tempList = []
            for repeatType in repeatClasses.allRepeats[repeatClass]: 
                wildcard='[N]*' + repeatType + '+'                                          # build wildcard string for N bases
                searchString = repeatType + minSize[repeatClass[:-7]] + wildcard            # concatenates mononuc[i] and minimum size for repeat type
                compiledRegEx = re.compile(searchString, re.IGNORECASE)
                tempList.append(compiledRegEx)
            if repeatClass[:-7]=='mononucleotide':
                self.mononucleotide = tempList
            elif repeatClass[:-7]=='dinucleotide':
                self.dinucleotide = tempList
            elif repeatClass[:-7]=='trinucleotide':
                self.trinucleotide = tempList
            elif repeatClass[:-7]=='tetranucleotide':
                self.tetranucleotide = tempList
            elif repeatClass[:-7]=='pentanucleotide':
                self.pentanucleotide = tempList
            elif repeatClass[:-7]=='hexanucleotide':
                self.hexanucleotide = tempList

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
    
    def genericMethod(self, passedSearchClass, repeatName, repeatUnits):
        """generic method for finding various microsatellite repeats"""
        # quick and dirty means of indexing position in passedSearchClass
        compiledRegExPos = 0                                                   
        sLength = len(self.seq)
        for compiledRegEx in passedSearchClass:
            iterator = compiledRegEx.finditer(self.seq)
            compIterator = compiledRegEx.finditer(self.revComplement(self.seq))
            name = repeatName + 'Letters'
            i = repeatClasses.allRepeats[name][compiledRegExPos]
            # list to help remove duplicates
            baseList = []
            # look for matches in forward sequence
            for match in iterator:
                bases = match.span()                                                
                # give start/end bases of repeat
                length = (bases[1] - bases[0]) / repeatUnits[repeatName]
                #self.msatResults[bases[0]+1] = ('%s repeat %s^%s found between bases %s to %s') % (repeatName.capitalize(), i, length, bases[0]+1, bases[1]+1)
                self.msatResults[bases[0]+1] = [bases[0]+1, i, length, bases[1]+1,repeatName.capitalize(), 'Forward']
                baseList.append(bases[0])
            # do the same on the reverse complement of the sequence
            for match in compIterator:                                              
                bases = match.span()
                # add if statement to remove repetitive finds...
                if sLength - bases[1] not in baseList:
                    length = (bases[1] - bases[0]) / repeatUnits[repeatName]
                    #seq = match.group()            
                    i = i.strip('()')
                    #self.msatResults[bases[0]+1] = ('Reverse complement of %s repeat %s, (%s)^%s found between bases %s to %s') % (repeatName, i, self.revComplement(i), length, sLength - bases[1]+1, sLength - bases[0]+1)
                    self.msatResults[bases[0]+1] = [sLength - bases[1]+1, (('(%s)') % (self.revComplement(i))), length, sLength - bases[0]+1, repeatName.capitalize(), (('Reverse:  complement of %s') % (i))]
            compiledRegExPos+=1

    def ephemeris(self, type, repeat):
        """Searches for microsatellite sequences (mononucleotide, dinucleotide, trinucleotide, tetranucleotide) in DNA string"""
        # we will store output for each repeat in dictionary keyed on start base #
        self.msatResults={}                                 
        # moved All (slow!) here to keep from double searching
        if 'All (slow!)' in type:                          
            self.genericMethod(repeat.mononucleotide, 'mononucleotide', repeat.units)
            self.genericMethod(repeat.dinucleotide,'dinucleotide', repeat.units)
            self.genericMethod(repeat.trinucleotide,'trinucleotide', repeat.units)
            self.genericMethod(repeat.tetranucleotide,'tetranucleotide', repeat.units)
            self.genericMethod(repeat.pentanucleotide,'pentanucleotide', repeat.units)
            self.genericMethod(repeat.hexanucleotide,'hexanucleotide', repeat.units)      
        else:
            for searchClass in type:
                if str(searchClass) == 'Mononucleotide':
                    self.genericMethod(repeat.mononucleotide, 'mononucleotide', repeat.units)        
                elif str(searchClass) == 'Dinucleotide':
                    self.genericMethod(repeat.dinucleotide,'dinucleotide', repeat.units)              
                elif str(searchClass) == 'Trinucleotide':
                    self.genericMethod(repeat.trinucleotide,'trinucleotide', repeat.units)
                elif str(searchClass) == 'Tetranucleotide':             
                    self.genericMethod(repeat.tetranucleotide,'tetranucleotide', repeat.units)
                elif str(searchClass) == 'Pentanucleotide':               
                    self.genericMethod(repeat.pentanucleotide,'pentanucleotide', repeat.units)
                elif str(searchClass) == 'Hexanucleotide':
                    self.genericMethod(repeat.hexanucleotide,'hexanucleotide', repeat.units)
        return self.msatResults
    
class excelSingleSpace:
    """class for csv module to work correctly"""
    delimiter = ','
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r'
#----------------------------------------------------------------------------------------
# program related functions for GUI (repeated within wx.Frame and other wx.Classes)
#----------------------------------------------------------------------------------------

class fileop:
    def __init__(self, outfile, infile):
        # get start time
        self.startTime = time.time()
        # open file for csv output
        self.file=open(outfile,'w')
        self.csvWriter = csv.writer(self.file, dialect = excelSingleSpace)
        # write some data
        self.csvWriter.writerow(['Clone','Start BP','Repeat','End BP','Type','Comments'])                                   
        # new-style biopython SeqIO iterator
        self.bioPythonHandle = open(infile,"rU")
 
    def printRunData(self, selection):
        for choice in selection:
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
        
    def beginSearch(self,selection,repeat,combine,primers):
        # initialize variable for number of sequence searched
        self.sequenceCount = 0
        # initialize list for sequence containing repeats
        self.overallRepeatSequences = []
        # overall repeats list
        self.overallRepeats = []
        for record in SeqIO.parse(self.bioPythonHandle,"fasta"):
            sequence = search(record.seq.tostring())
            dataOut=sequence.ephemeris(selection,repeat) 
            dictKeys=dataOut.keys()
            # sorts keys so bp locations will be in order by clone name
            dictKeys.sort()                                    
            if dictKeys:
                # add item for sequence containing repeat
                self.overallRepeatSequences.append(1)
                # writes dict values for sorted keys to output file
                for k in dictKeys:                                  
                    self.overallRepeats.append(1)
                    #dataList = dataOut[k].split()
                    #cloneResults = [record.id, ' '.join(dataList[:-7]), dataList[-7], ' '.join(dataList[-6:-3]), dataList[-3], dataList[-2], dataList[-1]]
                    cloneResults = [record.id, dataOut[k][0], (('%s^%s') % (dataOut[k][1], dataOut[k][2])), dataOut[k][3], dataOut[k][4], dataOut[k][5]]
                    self.csvWriter.writerow(cloneResults)
            elif not dictKeys:
                cloneResults = [record.id, "No repeats found"]        
                self.csvWriter.writerow(cloneResults)
            self.sequenceCount +=1       
        # print summary data IF NOT designing primers
        if not combine or not primers:
            self.printRunData(selection)
        # close open files
        self.bioPythonHandle.close()
        self.file.close()

    def genericError(self, errorItem):
        messageText = (('You must specify a %s') % (errorItem))
        windowTitle = 'Error'
        try:
            error = wx.MessageDialog(self.Parent, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        except:
            error = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        error.ShowModal()
        error.Destroy()
        
class prettify:
    def readMsatData(self, msatInfile):
        """get the freaking data"""
        dataDialect = csv.Sniffer().sniff(msatInfile)                        # try to determine dialect automagically
        #headerTrue = csv.Sniffer().has_header(msatInfile)                    # try to determine if file has header row
        msatData = open(msatInfile,'rU')                                               # open datafile for read
        msatReader = csv.reader(msatData, delimiter=',', dialect = dataDialect)     # setup csv.reader
        #if headerTrue:
        # skip header row
        msatReader.next()
        # use a list because indexes will not be unique (all the time...)
        # for clones w/ > 1 repeat array
        self.msatResults = []
        for row in msatReader:
            try:
                clone,start,repeat,end,reptype,comments = row
                # create as a tuple to allow sorting on > 1 key (1st then second (bp))
                self.msatResults.append((clone, int(start), repeat, int(end), reptype))
            except:
                pass
        self.msatResults.sort()
        msatData.close()
        
    def combine(self, infile, distance):
        #self.msatResults.sort()
        csvFile = open(infile,'w')
        combined_writer = csv.writer(csvFile, dialect = excelSingleSpace)
        combined_writer.writerow(['Clone','Start BP','Repeat','End BP','Type','Comments'])
        tempResults = []
        #pdb.set_trace()
        for item in range(len(self.msatResults)):
            # need to retain try:except as list is getting shorter
            try:
                i = 0
                primary = self.msatResults[item]
                while i <= len(self.msatResults):
                    # if last record, break
                    try:
                        secondary = self.msatResults[item+1]
                    except:
                        break
                    # if the names are the same... 
                    if primary[0] == secondary[0]:
                        name = 0
                        start = 1
                        repeat = 2
                        end = 3
                        # give compound repeat notation for differences of 0 BP
                        if primary[end] - secondary[start]==0:
                            newRepeat = primary[repeat]+secondary[repeat]
                            primary = (primary[name],primary[start],newRepeat,secondary[end], 'Compound')
                            # remove secondary
                            self.msatResults.remove(secondary)
                            i+=1
                        # give interrupted repeat notation for differences of <50 BP
                        elif primary[end] - secondary[start] < distance:
                            newRepeat = primary[repeat] + '...' + secondary[repeat]
                            primary = (primary[name], primary[start], newRepeat, secondary[end], 'Interrupted/Compound Interrupted')
                            # remove secondary
                            self.msatResults.remove(secondary)
                            i+=1
                    else:
                        break
                tempResults.append(primary)
                # write the record
                combined_writer.writerow(primary)
                #break
            except:pass
        csvFile.close()