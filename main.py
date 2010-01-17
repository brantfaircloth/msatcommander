#!/usr/bin/env python 

import os
import sys
import msat
import cPickle
import sqlite3
from Bio import SeqIO
from p3wrapr import primer
from PyQt4 import QtCore, QtGui
from ui_msatcommander import Ui_msatcommander

import pdb

class Window(QtGui.QWidget, Ui_msatcommander):
    '''stuff'''
    def __init__(self, parent = None):    
        QtGui.QWidget.__init__(self, parent)
        self.setupUi(self)
    
    def open(self):
        '''get the input filename - mutiple file input is not supported'''
        self.infile = QtGui.QFileDialog.getOpenFileName(self,\
            "Select Fasta to Search","~/", "Fasta (*.fsa *.fa *.fasta)")
        #TODO:Deal with other input formats e.g. fastq
        self.infileLength = open(self.infile, 'rU').read().count('>')
        return
        
    def save(self):
        '''get the directory to save the resutls'''
        self.outdir = QtGui.QFileDialog.getExistingDirectory(self,\
            "Select Directory to Save", "~/", QtGui.QFileDialog.ShowDirsOnly)
    
    def close(self):
        '''quit the program'''
        sys.exit()
    
    def accept(self):
        '''this is essentially the main loop'''
        self.open()
        self.save()
        self.createDbase()
        motifs = self.getMotifs()
        lengths = self.getLengths()
        self.generateCollection(motifs, lengths)
        self.readSearchSave()
        #TODO:  combined loci
        if self.designPrimersCheckBox.isChecked():
            self.designPrimers()
    
    def getMotifs(self):
        '''get the motifs for which we will search from user input'''
        return [m.isChecked() for m in (
            self.mononucCheckBox,
            self.dinucCheckBox,
            self.trinucCheckBox,
            self.tetranucCheckBox,
            self.pentanucCheckBox,
            self.hexanucCheckBox)
            ]
            
    def getLengths(self):
        '''get the lengths for the motifs for which we are searching'''
        return [m.value() for m in (
            self.mononucSpinBox,
            self.dinucSpinBox,
            self.trinucSpinBox,
            self.tetranucSpinBox,
            self.pentanucSpinBox,
            self.hexanucSpinBox)
            ]
    
    def generateCollection(self, motifs, lengths):
        '''generate the collection of regular expressions given the motifs 
        input and their requested lengths'''
        self.collection = ()
        possible_motifs = (
            msat.motif.mononucleotide,
            msat.motif.dinucleotide,
            msat.motif.trinucleotide,
            msat.motif.tetranucleotide,
            msat.motif.pentanucleotide,
            msat.motif.hexanucleotide
            )
        if self.perfectRepeatsCheckBox.isChecked():
            perfect = True
        else:
            perfect = False
        for m in enumerate(motifs):
            if m[1]:
                self.collection += (msat.seqsearch.MicrosatelliteMotif(possible_motifs[m[0]], lengths[m[0]], perfect),)
    
    def createDbase(self):
        '''create a database and associated tables to hold the results'''
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        self.conn = sqlite3.connect(os.path.join(str(self.outdir), 'msatcommander.sqlite'))
        self.cur = self.conn.cursor()
        # drop existing tables
        self.cur.execute('''DROP TABLE IF EXISTS records''')
        self.cur.execute('''DROP TABLE IF EXISTS sequences''')
        self.cur.execute('''DROP TABLE IF EXISTS microsatellites''')
        # create the sequence table
        self.cur.execute('''CREATE TABLE records (
        id int,
        name text,
        PRIMARY KEY(id)
        )''')
        self.cur.execute('''CREATE TABLE sequences (
        id int,
        seq blob,
        FOREIGN KEY(id) REFERENCES records(id)
        )''')
        # create the microsatellite table
        self.cur.execute('''CREATE TABLE microsatellites (
        id int,
        motif text,
        start int,
        end int,
        preceding int,
        following int,
        count int,
        FOREIGN KEY(id) REFERENCES records(id)
        )''')
        #TODO: move commit until after all operations?
        self.conn.commit()
    
    def searchForMotif(self, record, msat):
        '''Generalized microsatellite search function'''
        for repeat in range(len(msat.compiled)):
            temp_match = ()
            for m in msat.compiled[repeat].finditer(str(record.seq)):
                temp_match += ((m.span(),m.span()[0],len(record.seq)-m.span()[1]),)
            if temp_match:
                record.matches[msat.motif[repeat]] = temp_match
    
    def readSearchSave(self):
        '''read the infile, search the reads for msats'''
        index = 0
        # Setup progress bar - boo-yah!!
        self.pb = QtGui.QProgressDialog("Searching for microsatellites...",\
            "Cancel", 0, self.infileLength)
        self.pb.setWindowModality(QtCore.Qt.WindowModal)
        for record in SeqIO.parse(open(self.infile,'rU'), 'fasta'):
            # add matches attribute to record
            record.matches = {}
            # search for the motif in the record.seq; store results in matches
            # attribute, biatch
            for motif in self.collection:
                self.searchForMotif(record, motif)
            # pickle the sequence record, because we're gonna use it later
            # and insert it into sqlite (we have to run Binary on it, first)
            rPickle = cPickle.dumps(record, 1)
            self.cur.execute('''INSERT INTO records (id, name) VALUES (?,?)'''\
                , (index, record.name))
            self.cur.execute('''INSERT INTO sequences (id, seq) VALUES (?,?)'''\
            , (index, sqlite3.Binary(rPickle)))
            # go through the motifs and insert them to the dbase
            for match in record.matches:
                for motif in record.matches[match]:
                    count = (motif[0][1] - motif[0][0])/len(match)
                    self.cur.execute('''INSERT INTO microsatellites (id, \
                        motif, start, end, preceding, following, count) \
                        VALUES (?,?,?,?,?,?,?)''',(index, match, motif[0][0],\
                        motif[0][1], motif[1], motif[2], count))
            self.conn.commit()
            # we're manually indexing the primary key here
            index += 1
            # update the progress bar
            self.pb.setValue(index)
            # make sure that we abort if triggered
            if self.pb.wasCanceled():
                break
        self.pb.setValue(self.infileLength)
    
    def createPrimersTable(self):
        '''add tables to the dbase to hold the primers'''
        pass
    
    def designPrimers(self):
        '''design primers for those reads possessing msat repeats'''
        # setup basic primer design parameters
        settings = primer.Settings()
        settings.basic()
        # Update primer3 settings for mispriming library
        settings.params['PRIMER_MISPRIMING_LIBRARY'] = 'misprime_lib_weight'
        #TODO: override settings with user input
        
        
        # get the reads with msats from the dbase
        if not self.combineLociCheckBox.isChecked():
            self.cur.execute('''SELECT count(*) FROM microsatellites''')
            count = self.cur.fetchall()[0][0]
            self.cur.execute('''SELECT 
                sequences.id,
                sequences.seq,
                microsatellites.start,
                microsatellites.end
                FROM sequences, microsatellites
                WHERE sequences.id = microsatellites.id
                ''')
        else:
            self.cur.execute('''SELECT count(*) FROM combined_microsatellites''')
            count = self.cur.fetchall()[0][0]
            self.cur.execute('''SELECT sequences.id,
                sequences.seq,
                combined_microsatellites.start,
                combined_microsatellites.end
                FROM sequences, combined_microsatellites
                WHERE sequences.id = combined_microsatellites.id
                ''')
        sequences = self.cur.fetchall()
        for seq in sequences:
            # de-blob the sequence objects
            record = cPickle.loads(str(seq[1]))
            target = '%s,%s' % (seq[2], seq[3]-seq[2])
            primer3 = primer.Primers()
            primer3.pick(settings, sequence=str(record.seq), target=target, name = 'primers')
            try:
                print "PRIMER LEFT:  %s, PRIMER RIGHT: %s" % (primer3.primers[0]['PRIMER_LEFT_SEQUENCE'], primer3.primers[0]['PRIMER_LEFT_SEQUENCE'])
            except:
                pass
        QtCore.pyqtRemoveInputHook()
        pdb.set_trace()  


def qt_trace():
  '''Set a tracepoint in the Python debugger that works with Qt'''
  from PyQt4.QtCore import pyqtRemoveInputHook
  from pdb import set_trace
  pyqtRemoveInputHook()
  set_trace()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())