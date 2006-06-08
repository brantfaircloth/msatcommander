#!/usr/bin/python

#-----------------------------------------------------------------------------------------------|
# Microsatellite finder for Python:  an update to N. Dean Pentcheff's Ephemeris 1.0 for Perl    |
# available from http://www.uga.edu/srel/DNA_Lab/ephemeris%201.0.bin.                           |
#                                                                                               |
# Copyright (C) 2005 Brant C. Faircloth.  Modifications/recoding conducted solely by            |
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

class mods: 
    from string import translate
    """Class representing DNA as a string sequence.""" 
 
    def __init__(self, s): 
        """Create DNA instance initialized to string s.""" 
        self.seq = s 
    
    def complement(self):
        """return complementary dna sequence"""
        tab = string.maketrans('AGCTagct','TCGAtcga')
        output = string.translate(self.seq, tab)
        return output
    
class search: 
    from string import upper
    def __init__(self, s): 
        """Create DNA instance initialized to string s.""" 
        self.seq = s        
    def ephemeris(self):
        """Searches for microsatellite sequences (mononucleotide, dinucleotide, trinucleotide, tetranucleotide) in DNA string"""
        min_size={"mono":'{10,}',"di":'{7,}',"tri":'{5,}',"tetra":'{4,}'} #defines minimum size for repeat unit
        #defines repeat units (lowest alphabetical, unique, non-complementary) for which we are searching
        mononuc=['(A)','(C)'] 
        dinuc=['(AC)','(AG)','(AT)','(CG)'] 
        trinuc=['(AAC)','(AAG)','(AAT)','(ACC)','(ACG)','(ACT)','(AGC)','(AGG)','(ATC)','(CCG)']
        tetranuc=['(AAAC)','(AAAG)','(AAAT)','(AACC)','(AACG)','(AACT)','(AAGC)','(AAGG)','(AAGT)','(ACAG)','(ACAT)','(ACCC)','(ACCG)','(ACCT)','(ACGC)','(ACGT)','(ACTC)','(ACTG)','(AGAT)','(AATC)','(AATG)','(AATT)','(ACGG)','(AGCC)','(AGCG)','(AGGC)','(AGGG)','(ATCC)','(ATCG)','(ATGC)','(CCCG)','(CCGG)']         
        output={} #we will store output for each repeat type in a dictionary keyed on the starting base of the repeat
        for i in mononuc:
            wildcard='[N]*'+i+'+' #must build wildcard string for N bases
            search_string = i+min_size['mono']+wildcard #concatenates mononuc[i] and minimum size for repeat type
            test = re.compile(search_string, re.IGNORECASE) #compiles regex for concatenated values
            iterator = test.finditer(self.seq) #sets up iterative repeat finder
            iterator2 = test.finditer(mods(self.seq).complement()) #creates complement of sequence to search for complementary repeats
            for match in iterator:
                bases = match.span() #give the bases of the repeat
                length = bases[1]-bases[0] #gives length of repeat, multiplied by 10 here for mononucs
                output[bases[0]+1]=('Mononucleotide repeat %s^%s found between bases %s and %s.') % (i, length, bases[0]+1, bases[1]+1)
            for match in iterator2:
                bases = match.span()
                length = bases[1]-bases[0]
                seq=match.group()
                output[bases[0]+1]=('Reverse complement of mononucleotide repeat %s, (%s)^%s found between bases %s and %s.') % (i, string.upper(mods(seq[0]).complement()), length, bases[0]+1, bases[1]+1)      
        for i in dinuc:
            wildcard='[N]*'+i+'+' #must build wildcard string for N bases
            search_string = i+min_size['di']+wildcard
            test = re.compile(search_string, re.IGNORECASE)       
            iterator = test.finditer(self.seq)
            iterator2 = test.finditer(mods(self.seq).complement())
            for match in iterator:
                bases = match.span()
                length = ((bases[1]-bases[0])/2)
                output[bases[0]+1]=('Dinucleotide repeat %s^%s found between bases %s and %s.') % (i, length, bases[0]+1, bases[1]+1)
            for match in iterator2:
                bases = match.span()
                length = ((bases[1]-bases[0])/2)
                seq=match.group()
                output[bases[0]+1]=('Reverse complement of dinucleotide repeat %s, (%s)^%s found between bases %s and %s.') % (i, string.upper(mods(seq[0:2]).complement()), length, bases[0]+1, bases[1]+1)         
        for i in trinuc:
            wildcard='[N]*'+i+'+' #must build wildcard string for N bases
            search_string = i+min_size['tri']+wildcard
            test = re.compile(search_string, re.IGNORECASE)       
            iterator = test.finditer(self.seq)
            iterator2 = test.finditer(mods(self.seq).complement())
            for match in iterator:
                bases = match.span()
                length = ((bases[1]-bases[0])/3)
                output[bases[0]+1]=('Trinucleotide repeat %s^%s found between bases %s and %s.') % (i, length, bases[0]+1, bases[1]+1)
            for match in iterator2:
                bases = match.span()
                length = ((bases[1]-bases[0])/3)
                seq=match.group()
                output[bases[0]+1]=('Reverse complement of trinucleotide repeat %s, (%s)^%s found between bases %s and %s.') % (i, string.upper(mods(seq[0:3]).complement()), length, bases[0]+1, bases[1]+1)
        for i in tetranuc:
            wildcard='[N]*'+i+'+' #must build wildcard string for N bases
            search_string = i+min_size['tetra']+wildcard
            test = re.compile(search_string, re.IGNORECASE)       
            iterator = test.finditer(self.seq)
            iterator2 = test.finditer(mods(self.seq).complement())
            for match in iterator:
                bases = match.span()
                length = ((bases[1]-bases[0])/4)
                output[bases[0]+1]=('Tetranucleotide repeat %s^%s found between bases %s and %s.') % (i, length, bases[0]+1, bases[1]+1)
            for match in iterator2:
                bases = match.span()
                length = ((bases[1]-bases[0])/4)
                seq=match.group()
                output[bases[0]+1]=('Reverse complement of tetranucleotide repeat %s, (%s)^%s found between bases %s and %s.') % (i, string.upper(mods(seq[0:4]).complement()), length, bases[0]+1, bases[1]+1)
        return output

def getFiles(directory):
    fileList = [os.path.normcase(f) for f in os.listdir(directory)] #gets file name according to case sensitivity of file system
    fileList=[os.path.join(directory,f) for f in fileList if os.path.isfile(os.path.join(directory, f))]
                                                                #joins filenames with directory names for local path
    dsStore=os.path.join(directory, '.DS_Store')                #concats directory and .DS_Store
    if dsStore in fileList:
        fileList.remove(dsStore)                               #removes os x specific .DS_Store files
    return fileList                                            #returns file list to program

def getUserFiles():
    output=''
    optlist= getopt.getopt(sys.argv[1:], ':')
    input=optlist[1][1]
    input = os.path.abspath(string.strip(input))        #have to strip whitespace characters for dragging folders
    try:
        fileList=getFiles(input)                            #calls function above
    except:
        print 'No directory found, assuming single file input'
        fileList=input
        output=os.path.join(os.path.dirname(input),'output.txt')
    if not output:
        output=os.path.join(input,'output.txt')
    else:
        try:
            os.path.isdir(os.path.dirname(os.path.abspath(output)))
        except:
            print 'This is not a valid path.  Make sure you have entered the path correctly.'
            sys.exit()
    print ('\nOutput file written to %s') % (os.path.abspath(output))
    return fileList, output
    
def readInfo(files,output):
    file=open(output,'a') #opens file for output - append only to keep from overwriting
    file.write('Microsatellite repeats found in the following sequences: \n')
    if type(files) == str:       #for single file entries
        fileName=files
        f=open(fileName,'r') #opens files to read
        fileContents=f.read() #reads the bad boys
        lineEndings=['\r','\n'] #removes pesky line endings, if present
        for i in lineEndings:
            fileContents=fileContents.replace(i,'')
        dataOut=search(fileContents).ephemeris() #runs ephemeris method of search class to find SSRs
        dictKeys=dataOut.keys() #gets keys from dictionary returned from above
        dictKeys.sort() #sorts keys so bp locations will be in order
        file.write(('%s**********************%s**********************%s') % ('\n', fileName, '\n'))
        for k in dictKeys: #writes dict values for sorted keys to output file
            file.write(('%s %s') % (dataOut[k], '\n'))
    else:
        for i in files:
            fileName=i
            f=open(i,'r') #opens files to read
            fileContents=f.read() #reads the bad boys
            lineEndings=['\r','\n'] #removes pesky line endings, if present
            for i in lineEndings:
                fileContents=fileContents.replace(i,'')
            dataOut=search(fileContents).ephemeris() #runs ephemeris method of search class to find SSRs
            dictKeys=dataOut.keys() #gets keys from dictionary returned from above
            dictKeys.sort() #sorts keys so bp locations will be in order
            file.write(('%s**********************%s**********************%s') % ('\n', fileName, '\n'))
            for k in dictKeys: #writes dict values for sorted keys to output file
                file.write(('%s %s') % (dataOut[k], '\n'))
    file.close()

print '-----------------'
print (('If you would like to see more verbose output (what is going on) or have more control over input and output files and directories, use the command-line version of this program (microsat_command_line.py).  It is invoked with at least the following options %s%spython microsat_command_line.py -f /path/to/your/input/file %s') % ('\n','\n','\n'))
print '-----------------'
files=getUserFiles()
readInfo(files[0],files[1])