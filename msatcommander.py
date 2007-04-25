#!/usr/bin/env python
# encoding: utf-8
"""
msatcommander.py

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
import wx, mscfunc
import wx.xrc

GUI_FILENAME       = "mscGUI.xrc"
GUI_MAINFRAME_NAME = "mscGuiMainFrame"

class MyApp(wx.App):
    def OnInit(self):
        '''load controls'''
        self._do_layout()
        return True
    
    def _do_layout(self):
        '''read XRC file and populate variables'''
        self.res = wx.xrc.XmlResource(GUI_FILENAME)
        self.InitVar()
    
    def InitVar(self):
        '''initalize some other variables (non-GUI)'''
        self.selection              = [u"All (slow!)"]          # we will adjust this list later                             
        self.combineArraysFlag      = False                     # for safety (GUI default is false)
        self.designPrimersFlag      = False
        self.tagPrimersFlag         = False
        self.combDistance           = 50
        self.infile                 = False
        self.outfile                = False
        self.cag                    = True
        self.m13                    = True
        self.custom                 = None
        self.InitFrame()
    
    def InitFrame(self):
        '''initialize events and stuff'''
        self.frame = self.res.LoadFrame(None, GUI_MAINFRAME_NAME)
        self.mainPanel = wx.xrc.XRCCTRL(self.frame, 'mainPanel')
        # ============================== #
        # = msat search related values = #
        # ============================== #
        self.msatChoice = wx.xrc.XRCCTRL(self.mainPanel, "msatChoice")
        self.frame.Bind(wx.EVT_CHECKLISTBOX, self.onSelection, self.msatChoice)
        '''get length values for repeat classes (ugly)'''
        # mono length
        self.monoLength = wx.xrc.XRCCTRL(self.mainPanel, "monoLength")
        self.frame.Bind(wx.EVT_COMBOBOX, None, self.monoLength)
        #di length
        self.diLength = wx.xrc.XRCCTRL(self.mainPanel, "diLength")
        self.frame.Bind(wx.EVT_COMBOBOX, None, self.diLength)
        #tri length
        self.triLength = wx.xrc.XRCCTRL(self.mainPanel, "triLength")
        self.frame.Bind(wx.EVT_COMBOBOX, None, self.triLength)
        #tetra length
        self.tetraLength = wx.xrc.XRCCTRL(self.mainPanel, "tetraLength")
        self.frame.Bind(wx.EVT_COMBOBOX, None, self.tetraLength)
        #penta length
        self.pentaLength = wx.xrc.XRCCTRL(self.mainPanel, "pentaLength")
        self.frame.Bind(wx.EVT_COMBOBOX, None, self.pentaLength)
        #hexa length
        self.hexaLength = wx.xrc.XRCCTRL(self.mainPanel, "hexaLength")
        self.frame.Bind(wx.EVT_COMBOBOX, None, self.hexaLength)
        '''end length values for repeat classes'''        
        # ================================ #
        # = primer design related values = #
        # ================================ #
        self.frame.Bind(wx.EVT_CHECKBOX, self.onCombineArrays, id=wx.xrc.XRCID('combineArrays'))
        self.combDistanceBox = wx.xrc.XRCCTRL(self.mainPanel, "combinationDistance")
        self.frame.Bind(wx.EVT_TEXT, self.onCombineDistance, self.combDistanceBox)
        self.frame.Bind(wx.EVT_CHECKBOX, self.onDesignPrimers, id=wx.xrc.XRCID('designPrimers'))
        self.frame.Bind(wx.EVT_CHECKBOX, self.onTagPrimers, id=wx.xrc.XRCID('tagPrimers'))
        self.frame.Bind(wx.EVT_CHECKBOX, self.onCagTag, id=wx.xrc.XRCID('cagTag'))
        self.frame.Bind(wx.EVT_CHECKBOX, self.onM13Tag, id=wx.xrc.XRCID('m13Tag'))
        self.customBox = wx.xrc.XRCCTRL(self.mainPanel, 'customTag')
        self.frame.Bind(wx.EVT_TEXT, self.onCustom, self.customBox)
        # ============================================== #
        # = opening-, saving-, running-related events  = #
        # ============================================== #
        self.frame.Bind(wx.EVT_BUTTON, self.onOpen, id=wx.xrc.XRCID('fileOpen'))
        self.frame.Bind(wx.EVT_BUTTON, self.onSave, id=wx.xrc.XRCID('fileSave'))
        self.frame.Bind(wx.EVT_BUTTON, self.onRun, id=wx.xrc.XRCID('runProgram'))
        # frame display
        self.frame.Show()
    
    def onOpen(self,evt):
        messageText = 'Choose a File for Searching:'
        wildCard = "Text or Fasta files|*.txt;*.fsa;*.fasta|All files|*.*"
        dlg = wx.FileDialog(self.frame, messageText, wildcard=wildCard, style=wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.infile = dlg.GetPath()                                                       
            dlg.Destroy() 
        else:
            messageText = 'You must specify an input file!'
            windowTitle = 'Error'
            error = wx.MessageDialog(self.frame, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
            error.ShowModal()
            error.Destroy()        
            dlg.Destroy()
    
    def onSave(self, evt):
        if self.designPrimersFlag:
            # if designing primers, get directory for storage
            dlg = wx.DirDialog(self.frame, "Select A Directory:")
            if dlg.ShowModal() == wx.ID_OK:
                self.outDir = dlg.GetPath()
                fileName = '/microsatSearchOutput.csv'
                self.outfile = self.outDir + fileName
                dlg.Destroy()
        else:
            # else, get a file for output
            dlg = wx.FileDialog(self.frame, "Select A File:", style=wx.SAVE)
            if dlg.ShowModal() == wx.ID_OK:
                self.outfile = dlg.GetPath()
                dlg.Destroy()
                if self.outfile == self.infile:
                    messageText = 'You cannot overwrite the input file!'
                    windowTitle = 'Error'
                    self.warn(windowTitle, messageText)            
                    dlg.Destroy()
        
    def warn(self, windowTitle, messageText):
        error = wx.MessageDialog(self.frame, messageText, windowTitle, wx.OK | wx.ICON_INFORMATION)
        error.ShowModal()
        error.Destroy()
    
    def onRun(self,evt):
        '''run the program'''
        # check for infile and outfile
        if self.checkFiles() != False and self.checkSelection() != False:
            self.createPatterns()
            self.runSearch()
            if self.combineArraysFlag and not self.designPrimersFlag:
                self.combineArrays()
            if self.designPrimersFlag:
                self.getPrimers()
            if self.tagPrimersFlag:
                self.tagPrimers()
            else:
                self.runEnd()
    
    def combineArrays(self):
        '''combine msat arrays (for compound/interrupted repeats)'''
        messageText = 'Choose a location for the combined results file'
        windowTitle = 'Save combined array file as:'
        dlg = wx.FileDialog(self.frame, messageText, windowTitle, style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            self.prettyMsatOutfile = dlg.GetPath()
            dlg.Destroy()
            if self.prettyMsatOutfile == self.infile:
                messageText = 'You cannot overwrite the input file!'
                windowTitle = 'Error'
                self.warn(windowTitle,messageText)            
                dlg.Destroy()
        combined = mscfunc.prettify()
        combined.readMsatData(self.outfile)
        combined.combine(self.prettyMsatOutfile, self.combDistance)
    
    def tagPrimers(self):
        '''tag designed primers'''
        if not self.cag and not self.m13 and not self.custom:
                windowTitle = 'Error'
                messageText = 'You must choose a tagging option!'
                self.warn(windowTitle, messageText)
        import mscprimertag
        tagprimer = mscprimertag.primer3Tag(self.outDir,self.primerin)
        tagprimer.getData()
        #print 'custom tag is:', self.custom
        tagprimer.run(self.cag, self.m13, self.custom)
        self.runEnd()
    
    def getPrimers(self):
        '''design untagged primers'''
        import mscprimer
        primer = mscprimer.designPrimers(self.outDir)
        primer.readFastaData(self.infile)
        primer.readMsatData(self.outfile)
        if self.combineArraysFlag:
            primer.combine(self.combDistance)                   # create the combined file
        else:pass
        primer.header()
        primer.run()
        if self.tagPrimersFlag:
            # get name of untagged primer file
            self.primerin = primer.outfile
        primer.end()
    
    def checkFiles(self):
        if not self.infile or not self.outfile:
            windowTitle = "Error"
            messageText = "You have not selected an input/output file!"
            self.warn(windowTitle, messageText)
            return False
        else: return True
    
    def checkSelection(self):
        if not self.selection:
            windowTitle = 'Error'
            messageText = 'You have not selected any arrays'
            self.warn(windowTitle, messageText)
            return False
        else: return True
            
    def createPatterns(self):    
        # load the regex information (module imports repeatClasses)
        # create object
        self.repeat = mscfunc.repeats()
        # create search patterns
        self.repeat.createSearchPattern(
        self.monoLength.GetValue(),
        self.diLength.GetValue(),
        self.triLength.GetValue(),
        self.tetraLength.GetValue(),
        self.pentaLength.GetValue(),
        self.hexaLength.GetValue()
        )
        # patterns = repeat.class (i,e. repeat.mononucleotide)
    
    def runSearch(self):
        msat = mscfunc.fileop(self.outfile, self.infile)
        # send designPrimers to truncate summary data in outfile
        msat.beginSearch(self.selection, self.repeat, self.combineArraysFlag, self.designPrimersFlag)
        
    def runEnd(self):
        windowTitle = 'Search Completed'
        messageText = 'Finished scanning for microsatellite repeats'
        self.warn(windowTitle, messageText)
    
    def onSelection(self,evt):
        '''get selected values for microsat array search'''
        index = evt.GetSelection()
        label = self.msatChoice.GetString(index)
        if label not in self.selection:
            self.selection.append(label)
        elif label in self.selection:
            self.selection.remove(label)
    
    def onDesignPrimers(self,evt):
        '''set designPrimers = TRUE so we can create/save in a directory (not a file)'''
        index = evt.GetSelection()
        if index:
            self.designPrimersFlag = True
        else:
            self.designPrimersFlag = False
    
    def onTagPrimers(self,evt):
        index = evt.GetSelection()
        if index:
            self.tagPrimersFlag = True
        else:
            self.tagPrimersFlag = False
            
    def onCombineArrays(self,evt):
        index = evt.GetSelection()
        if index:
            self.combineArraysFlag = True
        else:
            self.combineArraysFlag = False
            
    def onCombineDistance(self,evt):
        try:
            self.combDistance = int(evt.GetString())
        except:
            windowTitle = 'Error'
            messageText = 'If you are combining repeats, you must give a integer distance (BP)'
            self.warn(windowTitle, messageText)
            
    def onCagTag(self,evt):
        index = evt.GetSelection()
        if index:
            self.cag = True
        else:
            self.cag = False
            
    def onM13Tag(self,evt):
        index = evt.GetSelection()
        if index:
            self.m13 = True
        else:
            self.m13 = False  
                    
    def onCustom(self,evt):
        self.custom = str(evt.GetString())
        if self.custom == '':
            self.custom = None
            self.customBox.SetValue('None')
            
if __name__ == '__main__':
    app = MyApp(False)
    app.MainLoop()
