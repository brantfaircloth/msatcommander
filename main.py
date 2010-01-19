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
    
    def message(self, text):
        self.message = QtGui.QMessageBox(self)
        self.message.setText(text)
        self.message.exec_()

    def checkPrimersOutput(self):
        '''check to make sure we've designed primers before we output them'''
        if not self.designPrimersCheckBox.isChecked():
            self.message = QtGui.QMessageBox.critical(self, 'Error', \
                '''You must (click) "Design Primers" to output primer \
sequences.''')
            self.primersCheckBox.setCheckState(QtCore.Qt.Unchecked)
        else:
            pass
            
    def checkTaggedPrimersOutput(self):
        '''check to make sure we've tagged primers before we output them'''
        if not self.tagPrimersCheckBox.isChecked() or not self.designPrimersCheckBox.isChecked():
            self.message = QtGui.QMessageBox.critical(self, 'Error',
                '''You must select "Design Primers" AND "Tag Primers" to \
output tagged primer sequences.''')
            self.taggedPrimersCheckBox.setCheckState(QtCore.Qt.Unchecked)
        else:
            pass
    
    def checkCombineRepeatsOutput(self):
        '''check to make sure we've combined loci before we output them'''
        if not self.combineLociCheckBox.isChecked():
            self.message = QtGui.QMessageBox.critical(self, 'Error',
                '''You must select "Combine Loci" to output combined \
repeats.''')
            self.combinedRepeatsCheckBox.setCheckState(QtCore.Qt.Unchecked)
        else:
            pass
    
    def checkRepeatsOutput(self):
        '''check to make sure we've searched for repeats before we combine'''
        if True not in [m.isChecked() for m in (
            self.mononucCheckBox,
            self.dinucCheckBox,
            self.trinucCheckBox,
            self.tetranucCheckBox,
            self.pentanucCheckBox,
            self.hexanucCheckBox)
            ]:
            self.message = QtGui.QMessageBox.critical(self, 'Error',
                '''You must select a repeat class for which to search before \
to output repeats.''')
            self.repeatsCheckBox.setCheckState(QtCore.Qt.Unchecked)
        else:
            pass
    
    def accept(self):
        '''this is essentially the main loop'''
        self.open()
        self.save()
        self.createDbase()
        motifs = self.getMotifs()
        lengths = self.getLengths()
        self.generateCollection(motifs, lengths)
        self.readSearchSave()
        #if self.combineLociCheckBox.isChecked()
        #    self.combineLoci()
        if self.designPrimersCheckBox.isChecked():
            self.designPrimers()
        self.outputResults()
    
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
        records_id int,
        id int,
        motif text,
        start int,
        end int,
        preceding int,
        following int,
        count int,
        FOREIGN KEY(records_id) REFERENCES records(id)
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
        msat_index = 0
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
                    self.cur.execute('''INSERT INTO microsatellites \
                        (records_id, id, motif, start, end, preceding, \
                        following, count) VALUES (?,?,?,?,?,?,?,?)''', \
                        (index, msat_index, match, motif[0][0],motif[0][1], \
                        motif[1], motif[2], count))
                    msat_index += 1
            self.conn.commit()
            # we're manually indexing the primary key here
            index += 1
            # update the progress bar
            self.pb.setValue(index)
            # make sure that we abort if triggered
            if self.pb.wasCanceled():
                break
        self.pb.setValue(self.infileLength)
    
    def createPrimersTable(self, name):
        '''add tables to the dbase to hold the primers'''
        # drop any current tables
        query = '''DROP TABLE IF EXISTS %s''' % name
        self.cur.execute(query)
        # create the new primers table
        query = ('''CREATE TABLE %s (
            records_id int,
            msats_id int,
            primer int,
            left text,
            left_sequence text,
            left_tm real,
            left_gc real,
            left_self_end real,
            left_self_any real,
            left_hairpin real,
            left_end_stability real,
            left_penalty real,
            right text,
            right_sequence text,
            right_tm real,
            right_gc real,
            right_self_end real,
            right_self_any real,
            right_hairpin real,
            right_end_stability real,
            right_penalty real,
            pair_product_size real,
            pair_compl_end real,
            pair_compl_any real,
            pair_penalty real,
            FOREIGN KEY(records_id) REFERENCES records(id),
            FOREIGN KEY(msats_id) REFERENCES microsatellites(id)
            )''' % name)
        self.cur.execute(query)
        self.conn.commit()
        
    def createTaggedPrimersTable(self):
        '''add tables to the dbase to hold the primers'''
        # drop any current tables
        self.cur.execute('''DROP TABLE IF EXISTS tagged_primers''')
        # create the new primers table
        self.cur.execute('''CREATE TABLE tagged_primers (
            records_id int,
            msats_id int,
            primer int,
            best int,
            tag text,
            tagged text,
            tag_seq test,
            common text,
            left text,
            left_sequence text,
            left_self_end real,
            left_self_any real,
            left_hairpin real,
            left_penalty real,
            right text,
            right_sequence text,
            right_self_end real,
            right_self_any real,
            right_hairpin real,
            right_penalty real,
            pair_product_size real,
            pair_compl_end real,
            pair_compl_any real,
            pair_penalty real,
            FOREIGN KEY(records_id) REFERENCES records(id),
            FOREIGN KEY(msats_id) REFERENCES microsatellites(id)
            )''')
        self.conn.commit() 
    
    def getMsatReads(self):
        '''get microsatellite containing reads from the database'''
        if not self.combineLociCheckBox.isChecked():
            self.cur.execute('''SELECT count(*) FROM microsatellites''')
            count = self.cur.fetchall()[0][0]
            self.cur.execute('''SELECT 
                sequences.id,
                sequences.seq,
                microsatellites.id,
                microsatellites.start,
                microsatellites.end
                FROM sequences, microsatellites
                WHERE sequences.id = microsatellites.records_id
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
        return count, sequences
    
    
    def storePrimers(self, table, record_id, msat_id, primers):
        '''store primers in the database'''
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        for i,p in primers.iteritems():
            if i != 'metadata':
                # create a copy of the dict, to which we add the
                # FOREIGN KEY reference
                td = p.copy()
                td['RECORD_ID'] = record_id
                td['MSAT_ID'] = msat_id
                td['PRIMER'] = i
                query = ('''INSERT INTO %s VALUES (
                    :RECORD_ID,
                    :MSAT_ID,
                    :PRIMER,
                    :PRIMER_LEFT,
                    :PRIMER_LEFT_SEQUENCE,
                    :PRIMER_LEFT_TM,
                    :PRIMER_LEFT_GC_PERCENT,
                    :PRIMER_LEFT_SELF_END_TH,
                    :PRIMER_LEFT_SELF_ANY_TH,
                    :PRIMER_LEFT_HAIRPIN_TH,
                    :PRIMER_LEFT_END_STABILITY,
                    :PRIMER_LEFT_PENALTY,
                    :PRIMER_RIGHT,
                    :PRIMER_RIGHT_SEQUENCE,
                    :PRIMER_RIGHT_TM,
                    :PRIMER_RIGHT_GC_PERCENT,
                    :PRIMER_RIGHT_SELF_END_TH,
                    :PRIMER_RIGHT_SELF_ANY_TH,
                    :PRIMER_RIGHT_HAIRPIN_TH,
                    :PRIMER_RIGHT_END_STABILITY,
                    :PRIMER_RIGHT_PENALTY,
                    :PRIMER_PAIR_PRODUCT_SIZE,
                    :PRIMER_PAIR_COMPL_END_TH,
                    :PRIMER_PAIR_COMPL_ANY_TH,
                    :PRIMER_PAIR_PENALTY
                    )''' % table)
                self.cur.execute(query, td)
        self.conn.commit()
    
    def storeTaggedPrimers(self, record_id, msat_id, primers, best=None):
        '''store primers in the database'''
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        for i,p in primers.iteritems():
            if i != 'metadata':
                # create a copy of the dict, to which we add the
                # FOREIGN KEY reference
                td = p.copy()
                td['RECORD_ID'] = record_id
                td['MSAT_ID'] = msat_id
                td['BEST'] = 0
                td['PRIMER'], td['TAG'], td['TAGGED'] = i.split('_')
                self.cur.execute('''INSERT INTO tagged_primers VALUES (
                    :RECORD_ID,
                    :MSAT_ID,
                    :PRIMER,
                    :BEST,
                    :TAG,
                    :TAGGED,
                    :PRIMER_TAG,
                    :PRIMER_TAG_COMMON_BASES,
                    :PRIMER_LEFT,
                    :PRIMER_LEFT_SEQUENCE,
                    :PRIMER_LEFT_SELF_END_TH,
                    :PRIMER_LEFT_SELF_ANY_TH,
                    :PRIMER_LEFT_HAIRPIN_TH,
                    :PRIMER_LEFT_PENALTY,
                    :PRIMER_RIGHT,
                    :PRIMER_RIGHT_SEQUENCE,
                    :PRIMER_RIGHT_SELF_END_TH,
                    :PRIMER_RIGHT_SELF_ANY_TH,
                    :PRIMER_RIGHT_HAIRPIN_TH,
                    :PRIMER_RIGHT_PENALTY,
                    :PRIMER_PAIR_PRODUCT_SIZE,
                    :PRIMER_PAIR_COMPL_END_TH,
                    :PRIMER_PAIR_COMPL_ANY_TH,
                    :PRIMER_PAIR_PENALTY
                    )''', (td))
        if best: # is the if necessary?
            best = best.split('_')
            self.cur.execute('''UPDATE tagged_primers 
                SET best = 1 WHERE 
                records_id = ? 
                AND msats_id = ?
                AND primer = ? 
                AND tag = ? 
                AND tagged = ?''',
                (record_id, msat_id, best[0], best[1], best[2]))
        self.conn.commit()    
    
    def designPrimers(self):
        '''design primers for those reads possessing msat repeats'''
        self.createPrimersTable('primers')
        # setup basic primer design parameters
        settings = primer.Settings()
        settings.basic()
        # Update primer3 settings to include the mispriming library
        settings.params['PRIMER_MISPRIMING_LIBRARY'] = 'misprime_lib_weight'
        #TODO: override settings with user input
        #
        #
        #
        if self.tagPrimersCheckBox.isChecked() or self.pigtailPrimersCheckBox.isChecked():
            self.createTaggedPrimersTable()
            if self.tagPrimersCheckBox.isChecked():
                # setup the settings for tagging primers
                tag_settings = primer.Settings()
                tag_settings.reduced(PRIMER_PICK_ANYWAY=1) 
        # get the reads with msats from the dbase
        count, sequences = self.getMsatReads()
        self.pb = QtGui.QProgressDialog("Searching for primers...",\
            "Cancel", 0, count)
        self.pb.setWindowModality(QtCore.Qt.WindowModal)
        index = 0
        for seq in sequences:
            # unpickle the sequence objects
            record = cPickle.loads(str(seq[1]))
            target = '%s,%s' % (seq[3], seq[4]-seq[3])
            primer3 = primer.Primers()
            primer3.pick(settings, sequence=str(record.seq), target=target, name = 'primers')
            # update dbase tables with primers
            if primer3.primers_designed:
                self.storePrimers('primers', seq[0], seq[2], primer3.primers)
                # if we're only pigtailing the designed primers
                if self.pigtailPrimersCheckBox.isChecked() and not self.tagPrimersCheckBox.isChecked():
                    primer3.pigtail(str(self.pigtailPrimersTagLineEdit.text()))
                    self.storeTaggedPrimers(seq[0], primer3.tagged_good, primer3.tagged_best_id)
            # if we are tagging
            if primer3.primers_designed and self.tagPrimersCheckBox.isChecked():
                cag, m13r, custom = None, None, None
                if self.cagTagCheckBox.isChecked(): cag = 'CAGTCGGGCGTCATCA'
                if self.m13rCheckBox.isChecked(): m13r = 'GGAAACAGCTATGACCAT'
                if self.customTagCheckBox.isChecked(): custom = str(self.customTagLineEdit.text())
                primer3.tag(tag_settings, CAG=cag, M13R=m13r, Custom = custom)
                if primer3.tagged_good:
                    if self.pigtailPrimersCheckBox.isChecked():
                        primer3.pigtail(str(self.pigtailPrimersTagLineEdit.text()))
                    #QtCore.pyqtRemoveInputHook()
                    #pdb.set_trace()
                    #update dbase tables with primers
                    self.storeTaggedPrimers(seq[0], seq[2], primer3.tagged_good, primer3.tagged_best_id)
            index += 1
            # update the progress bar
            self.pb.setValue(index)
            # make sure that we abort if triggered
            if self.pb.wasCanceled():
                break
        self.pb.setValue(count)
    
    
    def outputWriter(self, out, extension, header, data):
        f = open(out, 'w')
        if extension == 'csv':
            f.write('%s\n' % ','.join([str(x) for x in header])) 
            for d in data:
                f.write('%s\n' % ','.join([str(x) for x in d]))
        elif extension == 'tdt':
            f.write('%s\n' % ','.join([str(x) for x in header])) 
            for d in data:
                f.write('%s\n' % '\t'.join([str(x) for x in d]))
        f.close()
    
    def outputResults(self):
        '''write results to some sort of output file'''
        # set the output format
        if self.commaSeparatedRadioButton.isChecked():
            extension = 'csv'
        elif self.tabDelimitedRadioButton.isChecked():
            extension = 'tdt'
        
        if self.repeatsCheckBox.isChecked():
            out = 'msatcommander.microsatellites.%s' % extension
            self.cur.execute('SELECT * from microsatellites')
            data = self.cur.fetchall()
            self.cur.execute('PRAGMA table_info(microsatellites)')
            header = [x[1] for x in self.cur.fetchall()]
            self.outputWriter(out, extension, header, data)
        
        if self.combineLociCheckBox.isChecked():
            out = 'msatcommander.microsatellites.combined.%s' % extension
            self.cur.execute('SELECT * from combined_microsatellites')
        
        if self.primersCheckBox.isChecked() \
                and not self.tagPrimersCheckBox.isChecked() \
                and not self.pigtailPrimersCheckBox.isChecked():
            out = 'msatcommander.primers.%s' % extension
            # get the best unlabelled primers
            self.cur.execute('SELECT * from primers where primer = 0')
            data = self.cur.fetchall()
            self.cur.execute('PRAGMA table_info(primers)')
            header = [x[1] for x in self.cur.fetchall()]
            self.outputWriter(out, extension, header, data)
        
        if self.primersCheckBox.isChecked() \
                and self.pigtailPrimersCheckBox.isChecked() \
                and not self.tagPrimersCheckBox.isChecked():
            self.cur.execute('PRAGMA table_info(primers)')
            out = 'msatcommander.primers.%s' % extension
            header = [x[1] for x in self.cur.fetchall()]
            rows = ['primers.' + x for x in header]
            rows = ', '.join(rows)
            # get best primers based on best tagged primers
            query = '''SELECT %s FROM primers, tagged_primers WHERE 
                primers.records_id = tagged_primers.records_id 
                AND primers.msats_id = tagged_primers.msats_id
                AND primers.primer = tagged_primers.primer
                AND tagged_primers.best = 1''' % rows
            data = self.cur.execute(query)
            self.outputWriter(out, extension, header, data)
            self.cur.execute('PRAGMA table_info(tagged_primers)')
            out = 'msatcommander.pigtailed_primers.%s' % extension
            header = [x[1] for x in self.cur.fetchall()]
            self.cur.execute('SELECT * from tagged_primers where best = 1')
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
        
        if self.primersCheckBox.isChecked() \
                and self.taggedPrimersCheckBox.isChecked():
            self.cur.execute('PRAGMA table_info(primers)')
            out = 'msatcommander.primers.%s' % extension
            header = [x[1] for x in self.cur.fetchall()]
            rows = ['primers.' + x for x in header]
            rows = ', '.join(rows)
            # get best primers based on best tagged primers
            query = '''SELECT %s FROM primers, tagged_primers WHERE 
                primers.records_id = tagged_primers.records_id 
                AND primers.msats_id = tagged_primers.msats_id
                AND primers.primer = tagged_primers.primer
                AND tagged_primers.best = 1''' % rows
            data = self.cur.execute(query)
            self.outputWriter(out, extension, header, data)
            self.cur.execute('PRAGMA table_info(tagged_primers)')
            out = 'msatcommander.tagged_primers.%s' % extension
            header = [x[1] for x in self.cur.fetchall()]
            self.cur.execute('SELECT * from tagged_primers where best = 1')
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
        
        if self.taggedPrimersCheckBox.isChecked() \
                and not self.primersCheckBox.isChecked():
            self.cur.execute('PRAGMA table_info(tagged_primers)')
            out = 'msatcommander.tagged_primers.%s' % extension
            header = [x[1] for x in self.cur.fetchall()]
            self.cur.execute('SELECT * from tagged_primers where best = 1')
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())