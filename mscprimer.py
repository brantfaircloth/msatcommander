#!/usr/bin/env python
# encoding: utf-8
"""
mscprimer.py

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

import sys, os, csv, re, primer, subprocess, time
from Bio import SeqIO
#----------------------------------------------------------------------
# get files with repeats from excel spreadsheet output by msatCommander
#----------------------------------------------------------------------

class excelSingleSpace:
    """class for csv module to work correctly"""
    delimiter = ','
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r'

class primer3:
    """this class will run the functions associated with primer3_core"""
    
    def __init__(self, s): 
        """create DNA instance initialized to string s.""" 
        self.seq = s
        reg = re.compile('[^ACGTN]+')           # find non-ACGT bases
        self.seq = reg.sub('N',self.seq)        # convert non-ACGT bases to N
    
    def offset(self, start, length, offset):
        """function adds length to target region"""
        if start - offset <= 0: pass
        else: start -= offset                   # subtract 10 bases from target start (for non-specific crap)
        if length + offset > len(self.seq): pass
        else: length += offset*2                # add 20 bases to target length (for non-specific crap)
        self.start = start
        self.length = length                    
    
    def format(self, start, length, tempDir):
        """function builds temporary output file to be used by primer 3, subclassed below"""
        self.offset(start, length, 10)
        self.tempPrimer3File = tempDir + 'primer3_temp.txt'
        outfile=open(self.tempPrimer3File,'w')
        outfile.write(('SEQUENCE=%s\n') % (self.seq))
        outfile.write(('TARGET=%s,%s\n') % (self.start, self.length))
        features = primer.features
        for feat in features:
            outfile.write(('%s\n') % (feat))
        outfile.close()
        #------------------------
        # End primer3 file values
        #------------------------
    
    def search(self, clone, writer, outdir):
        primerInfo = outdir + clone
        primerInfo = re.sub('[|,.]+','_',primerInfo) + '_primer.txt'
        if os.name == 'posix':
            pathToPrettyPrimer3 = (('./primer3_core -format_output < \"%s\" > \"%s\"') % (self.tempPrimer3File, primerInfo))
            pathToPrimer3 = (('./primer3_core < \"%s\"') % (self.tempPrimer3File))
        elif os.name == 'nt':
            pathToPrettyPrimer3 = (('primer3_core -format_output < \"%s\" > \"%s\"') % (self.tempPrimer3File, primerInfo))
            pathToPrimer3 = (('primer3_core < \"%s\"') % (self.tempPrimer3File))
        # use new subprocess modules to allow things to run on windows - argh.
        #pdb.set_trace()
        subprocess.Popen(pathToPrettyPrimer3,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)    
        p=subprocess.Popen(pathToPrimer3,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True).stdout
        primer3Output = p.read()
        p.close()
        primer3Output = primer3Output.split('\n')
        #primer3Output=os.popen(pathToPrimer3).read().split('\n')
        primer0specs = {
        'PRIMER_RIGHT_SEQUENCE':None,
        'PRIMER_RIGHT_TM':None,
        'PRIMER_RIGHT_GC_PERCENT':None,
        'PRIMER_LEFT_SEQUENCE':None,
        'PRIMER_LEFT_TM':None,
        'PRIMER_LEFT_GC_PERCENT':None,
        'PRIMER_RIGHT_SELF_ANY':None,
        'PRIMER_RIGHT_SELF_END':None,
        'PRIMER_RIGHT_END_STABILITY':None,
        'PRIMER_RIGHT_PENALTY':None,
        'PRIMER_RIGHT_EXPLAIN':None,
        'PRIMER_LEFT_SELF_ANY':None,
        'PRIMER_LEFT_SELF_END':None,
        'PRIMER_LEFT_END_STABILITY':None,
        'PRIMER_LEFT_PENALTY':None,	
        'PRIMER_LEFT_EXPLAIN':None,
        'PRIMER_LEFT':None,
        'PRIMER_RIGHT':None
        }
        primer1specs,primer2specs,primer3specs,primer4specs = primer0specs.copy(),primer0specs.copy(),primer0specs.copy(),primer0specs.copy()
        listOfDicts = [primer0specs,primer1specs,primer2specs,primer3specs,primer4specs]
        compileStatements = 'PRIMER_(RIGHT|LEFT)[A-Z1-9_]*='
        i=0
        for primerDict in listOfDicts:
            #reg = re.compile(compileStatements[i])
            reg = re.compile(compileStatements)
            #print 'primer 3 output', primer3Output
            for item in primer3Output:
                m = reg.match(item)
                if m:
                    #m.group()
                    splitItem = item.split('=')
                    primerDict[splitItem[0]]=splitItem[1]
                else: pass
            i+=1
        # get primer start positions
        try:
            upperStart = primer0specs['PRIMER_LEFT'].split(',')[0]
            lowerStart = primer0specs['PRIMER_RIGHT'].split(',')[0]            
        except:
            upperStart = None
            lowerStart = None
        outputdata = [
            clone,
            primer0specs['PRIMER_LEFT_SEQUENCE'],
            primer0specs['PRIMER_LEFT_TM'],
            primer0specs['PRIMER_LEFT_GC_PERCENT'],
            upperStart,
            primer0specs['PRIMER_RIGHT_SEQUENCE'],
            primer0specs['PRIMER_RIGHT_TM'],
            primer0specs['PRIMER_RIGHT_GC_PERCENT'],            
            lowerStart
            ]
        if outputdata[1]:
            writer.writerow(outputdata)
        else:
            pass
        # remove tempfile
        try:
            os.remove(self.tempPrimer3File)
        except OSError:
            time.sleep(0.1)
            os.remove(self.tempPrimer3File)
        
class designPrimers:
    def __init__(self, outDir):
        self.msatDir = outDir
        if os.name == 'posix':
            self.outDir = outDir + '/primer3/'
        elif os.name == 'nt':
            self.outDir = outDir + '\\primer3\\'
        try:
            os.mkdir(self.outDir)
        except:
            pass
        self.outfile = self.outDir + 'primerNoTag.csv'
        self.csvFile = open(self.outfile,'w')
        self.outfile_writer = csv.writer(self.csvFile, dialect = excelSingleSpace)

    def readFastaData(self, fastaInfile):
        bioPythonHandle = open(fastaInfile,"rU")
        self.fastaDict = {}
        for record in SeqIO.parse(bioPythonHandle,"fasta"):
            sequence = record.seq.tostring()
            title = record.id
            self.fastaDict[title]=sequence
        bioPythonHandle.close()
    
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
        
    def combine(self, distance):
        #self.msatResults.sort()
        tagOutfile = self.msatDir + '/primerCombined.csv'
        csvFile = open(tagOutfile,'w')
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
        #replace self.msatResults with tempResults
        self.msatResults = tempResults
    
    def header(self):
        headerRow = [
            'CLONE',
            'PRIMER_LEFT',
            'LEFT_TM',
            'LEFT_GC_PERCENT',
            'LEFT_START',
            'PRIMER_RIGHT',
            'RIGHT_TM',
            'RIGHT_GC_PERCENT',
            'RIGHT_START'
            ]
        self.outfile_writer.writerow(headerRow)
    
    def run(self):
        for clone in self.msatResults:
            #print searchResults[clone]
            name = clone[0]
            length = clone[3] - clone[1]
            start = clone[1]
            seq = primer3(self.fastaDict[name])
            seq.format(start, length, self.outDir)
            seq.search(name, self.outfile_writer, self.outDir)
            #os.remove('primer3_temp.txt')
    
    def end(self):
        self.csvFile.close()
