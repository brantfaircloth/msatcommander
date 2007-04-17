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

import string, re, wx, time, csv, os.path, sys
#from Bio import Fasta
from Bio.SeqIO import SequenceIterator
 
idAbout = wx.NewId()
idOpen  = wx.NewId()
idSave  = wx.NewId()
idQuit  = wx.NewId()
idSearch = wx.NewId()
idConversion = wx.NewId()
idExit  = wx.NewId()
idListBox = wx.NewId()

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
    
    def genericMethod(self,passedSearchClass, repeat):
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
        
        """Searches for microsatellite sequences (mononucleotide, dinucleotide, trinucleotide, tetranucleotide) in DNA string"""        

        # we will store output for each repeat in dictionary keyed on start base #
        self.msatResults={}                                 
        # moved All (slow!) here to keep from double searching
        if 'All (slow)' in type:                          
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
#----------------------------------------------------------------------------------------
# program related functions for GUI (repeated within wx.Frame and other wx.Classes)
#----------------------------------------------------------------------------------------

class fileFunctions:
    
    def __init__(self, Parent):
        """ get-set parent window """
        self.Parent = Parent  
    
    def OnOpen(self, event):
        messageText = 'Choose a File for Searching'
        wildCard = "Text or Fasta files|*.txt;*.fsa;*.fasta|All files|*.*"
        try:
            dlg = wx.FileDialog(self.Parent, messageText, wildcard=wildCard)
        except:
            dlg = wx.FileDialog(self, messageText, wildcard=wildCard)
        dlg.SetStyle(wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.infile = dlg.GetPath()                                                       
            #self.userChoice()
            #print 'infile is...', self.infile
            dlg.Destroy()
            #return infile  
        else:
            messageText = 'You must specify an input file!'
            windowTitle = 'Error'
            try:
                error = wx.MessageDialog(self.Parent, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
            except:
                error = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
            error.ShowModal()
            error.Destroy()        
            dlg.Destroy()
            #self.OnOpen(event=None)  # for extra repetitive window opening - e.g. until input

    def SaveFile(self, event):
        messageText = 'Choose a Location to Save the Output File'
        wildCard = "Comma-Separated|*.csv"
        defDirectory = os.path.dirname(self.infile)
        try:
            dlg = wx.FileDialog(self.Parent, messageText, wildcard=wildCard, defaultDir=defDirectory)
        except:
            dlg = wx.FileDialog(self, messageText, wildcard=wildCard, defaultDir=defDirectory)
        dlg.SetStyle(wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.outfile = dlg.GetPath()
            dlg.Destroy()
            if self.outfile == self.infile:
                messageText = 'You cannot overwrite the input file!'
                windowTitle = 'Error'
                try:
                    error = wx.MessageDialog(self.Parent, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
                except:
                    error = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
                error.ShowModal()
                error.Destroy()            
                dlg.Destroy()
                self.SaveFile(event=None)  # for extra repetitive window opening - e.g. until input
            #return outfile
        else:
            messageText = 'You must specify an output file!'
            windowTitle = 'Error'
            try:
                error = wx.MessageDialog(self.Parent, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
            except:
                error = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
            error.ShowModal()
            error.Destroy()            
            dlg.Destroy()  
            #self.SaveFile(event=None)  # for extra repetitive window opening - e.g. until input
            
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
    
    def runEnd(self):
        messageText = 'Finished scanning for microsatellite repeats'
        windowTitle = 'Search Completed'
        dlg = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
    
    def RunSearch(self, event):
        self.startTime = time.time()
        # opens file for output (overwrite existing)
        file=open(self.outfile,'w')                                   
        self.csvWriter = csv.writer(file, dialect = excelSingleSpace)
        self.csvWriter.writerow(['Clone','Repeat Info','Repeat Count','Location','Start BP','','End BP'])
        # new-style biopython SeqIO iterator
        handle = open(self.infile,"rU")
        # initialize variable for number of sequence searched
        self.sequenceCount = 0
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
        self.printRunData()
        # close open files
        handle.close()
        file.close()
        self.runEnd()

    def genericError(self, errorItem):
        messageText = (('You must specify a %s') % (errorItem))
        windowTitle = 'Error'
        try:
            error = wx.MessageDialog(self.Parent, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        except:
            error = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        error.ShowModal()
        error.Destroy()
        
class fileConversions:
    
    def RunConversion():
        pass

#---------------------------------
# Main class for wxPython interface
#---------------------------------
class msatCommand(wx.App):
    def OnInit(self):
        frame = msatCommandFrame(None)
        frame.Show(True)
        self.SetTopWindow(frame)
        return True

class DemoPanel(wx.Panel, fileFunctions):
    """Panel gives main user interface for program - parent is frame"""
    def __init__(self, Parent, *args, **kwargs):
        """Create the DemoPanel."""
        wx.Panel.__init__(self, Parent, *args, **kwargs)
        # get parent of window (frame class)
        self.Parent = Parent
        # label the panel
        wx.StaticText(Parent, -1, pos=(5,5), label='Choose Repeat Class(es):')
        
        #present user with a open file button in GUI
        fileSelect = wx.Button(Parent, label="Select file to scan...", size=(200, 20), pos=(220,40))
        fileSelect.Bind(wx.EVT_BUTTON, self.OnOpen)       
        #present user with a save file button in GUI
        fileSave =  wx.Button(Parent, label="Select file for output...", size=(200, 20), pos=(220,90))
        fileSave.Bind(wx.EVT_BUTTON, self.SaveFile )
        
        # Create a checklistbox on an existing panel widget
        self.checklist = wx.CheckListBox(Parent, pos=(5, 25), size=(200,150))
        # Present user with a list of choices in the GUI (these begin w/ capital here but not in variable names)
        list = ['Mononucleotide',
        'Dinucleotide',
        'Trinucleotide',
        'Tetranucleotide',
        'Pentanucleotide', 
        'Hexanucleotide', 
        'All (slow)'
        ]
        # fire the items into the checklistbox
        self.checklist.InsertItems(items=list, pos=0)
        
        # Bind checklist to event
        self.checklist.Bind(wx.EVT_CHECKLISTBOX, self.OnSelection)
        
        # present user with a run search button in GUI
        runProgram =  wx.Button(Parent, label="Search file", size=(200, 20), pos=(220,140))
        runProgram.Bind(wx.EVT_BUTTON, self.RunSearch)
        
        # label options area
        wx.StaticText(Parent, -1, pos=(5,180), label='Options:')
        # add a checkbox to allow user to show clone without repeats
        noRepeats = wx.CheckBox(Parent, label="Show clones with no repeats", pos=(14,200))
        self.noRepeats = False
        noRepeats.Bind(wx.EVT_CHECKBOX, self.noRepeatFunction)
        
    def OnSelection(self,event):
        checked = []
        for i in range(self.checklist.GetCount()):
            if self.checklist.IsChecked(i):
                checked.append(self.checklist.GetString(i))       
        self.selection = checked
        return self.selection
    
    def noRepeatFunction(self,event):
        self.noRepeats = True
        return self.noRepeats

          
        
class msatCommandFrame(wx.Frame, fileFunctions):
    """Main wxPython frame for msatCommand"""
    title = "msatCommand"
    def __init__(self, parent):
        wx.Frame.__init__(self,parent,-1,self.title,size=(450,285), style=wx.DEFAULT_FRAME_STYLE)
        color=(255,255,255)
        self.SetBackgroundColour(color)
        self.CreateStatusBar()
        self.SetStatusText("")
        
        #create the file menu
        menu1 = wx.Menu()
        menu1.Append(idOpen, "&Open...\tCtrl-O","Open a file for searching")
        menu1.Append(idSave, "&Save output as...\tCtrl-S","Save output file as...")
        menu1.Append(idSearch, "&Search File\tCtrl-Z", "Convert open file")
        menu1.Append(idQuit, "&Quit\tCtrl-Q","Quit")

        #and the conversion menu
        #menu2 = wx.Menu()
       # menu2.Append(idConversion, "&Fasta Conversion\tCtrl-F", "Convert input file format")
        
        #and the help menu
        menu3 = wx.Menu()
        menu3.Append(idAbout, "&About\tCtrl-A", "Display information about this program")
        
        #add menu items to menubar
        menuBar = wx.MenuBar()
        menuBar.Append(menu1, "&File")
        #menuBar.Append(menu2, "&Conversion")
        menuBar.Append(menu3, "&Help") 
        self.SetMenuBar(menuBar)
        
        #create a panel
        self.Panel = DemoPanel(self)
        
        #setup events
        wx.EVT_MENU(self, idOpen, self.OnOpen)
        wx.EVT_MENU(self, idSave, self.SaveFile)
        wx.EVT_MENU(self, idSearch, self.SearchFunc)
        wx.EVT_MENU(self, idAbout, self.OnAbout)
        wx.EVT_MENU(self, idQuit, self.Quit)
        
        #self.Bind(wx.EVT_MENU, self.OnOpen)
        #self.Bind(wx.EVT_MENU, self.SaveFile)
        #self.Bind(wx.EVT_MENU, self.RunSearch, id=idSearch)
        #self.Bind(wx.EVT_MENU, fileConversions.RunConversion, id=idConversion)
        #self.Bind(wx.EVT_MENU, self.OnAbout, id=idAbout)        
        #self.Bind(wx.EVT_MENU, self.OnExit, id=idExit)
        
    
    def SearchFunc(self,event):
        #pass
        DemoPanel(self).OnSelection
        #print 'got selection'
        self.selection = DemoPanel(self).OnSelection
        #print self.selection
        self.RunSearch
    
    def OnAbout(self, event):
        messageText = 'This program searches Fasta files for microsatellite repeats'
        windowTitle = 'msatCommand'
        dlg = wx.MessageDialog(self, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        
       
    def OnExit(self, event):
        self.Close()
        #self.Close(true)
    
    def Quit(self,event):
        sys.exit()
        


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
    app = msatCommand(redirect=True)
    app.MainLoop()