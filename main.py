#!/usr/bin/env python 

import os
import sys
import cPickle
from pysqlite2 import dbapi2 as sqlite3
import operator
#import multiprocessing
from Bio import SeqIO
import msat
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
            "Select Fasta to Search",os.path.expanduser('~'), "Fasta (*.fsa *.fa *.fasta)")
        # ==================
        # = DEBUG EASIFIER =
        # ==================
        #self.infile = QtGui.QFileDialog.getOpenFileName(self,\
        #    "Select Fasta to Search",os.getcwd(), "Fasta (*.fsa *.fa *.fasta)")
        #TODO:Deal with other input formats e.g. fastq
        if self.infile:
            self.infileLength = open(self.infile, 'rU').read().count('>')
        else:
            self.infile = None
        return
        
    def save(self):
        '''get the directory to save the resutls'''
        self.outdir = QtGui.QFileDialog.getExistingDirectory(self,\
            "Select Directory to Save", "~/", QtGui.QFileDialog.ShowDirsOnly)
        if not self.outdir:
            self.outdir = None
    
    def close(self):
        '''quit the program'''
        sys.exit()
    
    def message(self, text):
        self.message = QtGui.QMessageBox(self)
        self.message.setText(text)
        self.message.exec_()

    
    def clickCombineLociCheckBox(self):
        if self.combineLociCheckBox.isChecked():
            if not self.combinedRepeatsCheckBox.isChecked():
                self.combinedRepeatsCheckBox.setCheckState(QtCore.Qt.Checked)
        else:
            self.combinedRepeatsCheckBox.setCheckState(QtCore.Qt.Unchecked)
        
    
    def clickRepeatsCheckBox(self):
        '''when you click a motif check box, set the state of the repeatsCheckBox'''
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        if self.sender().isChecked():
            if not self.repeatsCheckBox.checkState():
                self.repeatsCheckBox.setCheckState(QtCore.Qt.Checked)
        else:
            if True not in [m.isChecked() for m in (
                self.mononucCheckBox,
                self.dinucCheckBox,
                self.trinucCheckBox,
                self.tetranucCheckBox,
                self.pentanucCheckBox,
                self.hexanucCheckBox)
                ]:
                self.repeatsCheckBox.setCheckState(QtCore.Qt.Unchecked)
    
    def clickDesignPrimersCheckBox(self):
        '''when you click the design primers check box, select to ouput primer files'''
        if self.designPrimersCheckBox.isChecked():
            if not self.primersCheckBox.checkState():
                self.primersCheckBox.setCheckState(QtCore.Qt.Checked)
        else:
            self.primersCheckBox.setCheckState(QtCore.Qt.Unchecked)
    
    def clickTagPrimersCheckBox(self):
        '''when you click the design primers check box, select to output 
        primer files.  When unclicked, deselect what was selected.'''
        # if it's being checked
        if self.tagPrimersCheckBox.isChecked():
            if not self.designPrimersCheckBox.checkState():
                self.designPrimersCheckBox.setCheckState(QtCore.Qt.Checked)
                self.primersCheckBox.setCheckState(QtCore.Qt.Checked)
            if not self.taggedPrimersCheckBox.checkState():
                self.taggedPrimersCheckBox.setCheckState(QtCore.Qt.Checked)
                self.cagTagCheckBox.setCheckState(QtCore.Qt.Checked)
                self.m13rCheckBox.setCheckState(QtCore.Qt.Checked)
        else:
            self.taggedPrimersCheckBox.setCheckState(QtCore.Qt.Unchecked)
            self.cagTagCheckBox.setCheckState(QtCore.Qt.Unchecked)
            self.m13rCheckBox.setCheckState(QtCore.Qt.Unchecked)
    
    def clickKeepDatabaseCheckBox(self):
        '''when you unclick the keepDatabaseCheckBox, unclick the 
        keepSequenceRecords check box'''
        if self.keepDatabaseCheckBox.isChecked():
            pass
        else:
            self.keepDbaseSequenceRecords.setCheckState(QtCore.Qt.Unchecked)
    
    def checkPrimersOutput(self):
        '''check to make sure we've designed primers before we output them'''
        if not self.designPrimersCheckBox.isChecked():
            self.message = QtGui.QMessageBox.critical(self, 'Error', \
                '''You must (click) "Design Primers" to output primer \
sequences.''')
            self.primersCheckBox.setCheckState(QtCore.Qt.Unchecked)
        else:
            pass
            
    def checkKeepDatabaseCheckBox(self):
        '''check to make sure if we are keeping sequence reads, that we've
        also decided to keep the dbase'''
        if not self.keepDatabaseCheckBox.isChecked():
            self.message = QtGui.QMessageBox.critical(self, 'Error',
                '''You must keep the database before you can keep the sequence \
records.''')
            self.keepDbaseSequenceRecords.setCheckState(QtCore.Qt.Unchecked)
    
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
        if True not in [m.isChecked() for m in (
            self.mononucCheckBox,
            self.dinucCheckBox,
            self.trinucCheckBox,
            self.tetranucCheckBox,
            self.pentanucCheckBox,
            self.hexanucCheckBox)
            ]:
            self.message = QtGui.QMessageBox.critical(self, 'Error',
                '''You have not chosen to search for anything''')
        else:
            self.open()
            if self.infile:
                self.save()
                if self.outdir:
                    self.createDbase()
                    motifs = self.getMotifs()
                    lengths = self.getLengths()
                    self.generateCollection(motifs, lengths)
                    self.readSearchSave()
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
        if self.keepDatabaseCheckBox.isChecked():
            self.conn = sqlite3.connect(os.path.join(str(self.outdir), 'msatcommander.sqlite'))
        else:
            self.conn = sqlite3.connect(":memory:")
        self.cur = self.conn.cursor()
        # drop existing tables
        self.cur.execute('''DROP TABLE IF EXISTS combined_components''')
        self.cur.execute('''DROP TABLE IF EXISTS sequences''')
        self.cur.execute('''DROP TABLE IF EXISTS microsatellites''')
        self.cur.execute('''DROP TABLE IF EXISTS primers''')
        self.cur.execute('''DROP TABLE IF EXISTS tagged''')
        self.cur.execute('''DROP TABLE IF EXISTS combined''')
        self.cur.execute('''DROP TABLE IF EXISTS records''')
        # turn on foreign key support
        self.cur.execute('''PRAGMA foreign_keys = ON''')
        # create the sequence table
        self.cur.execute('''CREATE TABLE records (
        id int,
        name text,
        PRIMARY KEY(id)
        )''')
        if self.keepDbaseSequenceRecords.isChecked():
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
        if self.combineLociCheckBox.isChecked():
            self.cur.execute('''CREATE TABLE combined (
            records_id int,
            id int,
            motif text,
            start int,
            end int,
            preceding int,
            following int,
            members int,
            PRIMARY KEY(records_id, id)
            FOREIGN KEY(records_id) REFERENCES records(id)
            )''')
            self.cur.execute('''CREATE TABLE combined_components (
            records_id int,
            combined_id int,
            motif text,
            length int,
            FOREIGN KEY(records_id, combined_id) REFERENCES combined(records_id, id)
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
    
    def addLocus(self, combined, bt, ct):
        '''convenience function to add loci to be combined to the combined tuple'''
        if not bt in combined:
            combined += (bt,)
        if not ct in combined:
            combined += (ct,)
        return combined
    
    def combineLoci(self, record, min_distance):
        '''combined adjacent loci - this is somewhat cumbersome due to the
        format of the matches returned from msat (a dict with keys = motif).
        Essentially, we are running a pairwise comparison across all motifs
        located to determine which are within a predetermined distance from 
        one another'''
        temp_combined = []
        reorder = ()
        # turn our dict into something more useful for this purpose
        for motif in record.matches:
            for pos,val in enumerate(record.matches[motif]):
                reorder += ((motif, pos, val[0][0], val[0][1], val[1], val[2]),)
        # sort it
        reorder = sorted(reorder, key=operator.itemgetter(2))
        # combine adjacent loci at < min_distance
        for i in reorder:
            included = False
            if not temp_combined:
                temp_combined.append([i])
            else:
                for gp, g in enumerate(temp_combined):
                    for elem in g:
                        if i[2] - elem[3] <= min_distance:
                            temp_combined[gp].append(i)
                            included = True
                            break
                if not included:
                    temp_combined.append([i])
        # re-key
        for group in temp_combined:
            motifs = []
            if len(group) > 1:
                gs = group[0][2]
                ge = group[-1][2]
                gp = group[0][4]
                gf = group[-1][5]
            else:
                gs, ge = group[0][2], group[0][3]
                gp, gf = group[0][4], group[0][5]
            name = ''
            member_count = 0
            for pos,member in enumerate(group):
                if pos + 1 < len(group):
                    dist = group[pos + 1][3] - group[pos][3]
                    if dist > 1:
                        spacer = '...'
                    else:
                        spacer = ''
                else:
                    spacer = ''
                length = (member[3]-member[2])/len(member[0])
                name += '%s(%s)%s' % (member[0], length, spacer)
                motifs.append([member[0],length])
                member_count += 1
            record.combined[name] = (((gs, ge), gp, gf, member_count, motifs),)
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
        return record
    
    def readSearchSave(self):
        '''read the infile, search the reads for msats'''
        # setup multiprocessing
        index           = 0
        msat_index      = 0
        combine_index   = 0
        # Setup progress bar - boo-yah!!
        self.pb = QtGui.QProgressDialog("Searching for microsatellites...",\
            "Cancel", 0, self.infileLength)
        self.pb.setWindowModality(QtCore.Qt.WindowModal)
        if self.designPrimersCheckBox.isChecked():
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            # setup basic primer design parameters
            settings = primer.Settings()
            settings.basic(path=os.path.join(os.getcwd(), 'primer3_config/'))
            # Update primer3 settings to include the mispriming library
            settings.params['PRIMER_MISPRIMING_LIBRARY'] = 'misprime_lib_weight'
            # Update the primer3 settings with user choices/defaults:
            settings.params['PRIMER_PRODUCT_SIZE_RANGE'] = str(self.primerProductSizeTextBox.text())
            settings.params['PRIMER_MIN_TM']             = float(self.primerMinTmSpinBox.value())
            settings.params['PRIMER_OPT_TM']             = float(self.primerOptTmSpinBox.value())
            settings.params['PRIMER_MAX_TM']             = float(self.primerMaxTmSpinBox.value())
            settings.params['PRIMER_MIN_SIZE']           = int(self.primerMinSizeSpinBox.value())
            settings.params['PRIMER_OPT_SIZE']           = int(self.primerOptSizeSpinBox.value())
            settings.params['PRIMER_MAX_SIZE']           = int(self.primerMaxSizeSpinBox.value())
            settings.params['PRIMER_MIN_GC']             = float(self.primerMinGcSpinBox.value())
            settings.params['PRIMER_MAX_GC']             = float(self.primerMaxGcSpinBox.value())
            settings.params['PRIMER_MAX_POLY_X']         = int(self.primerMaxPolyXSpinBox.value())
            settings.params['PRIMER_MAX_SELF_ANY_TH']           = float(self.primerMaxSelfAnySpinBox.value())
            settings.params['PRIMER_PAIR_MAX_COMPL_ANY_TH']     = float(self.primerMaxPairAnySpinBox.value())
            settings.params['PRIMER_MAX_SELF_END_TH']           = float(self.primerMaxSelfEndSpinBox.value())
            settings.params['PRIMER_PAIR_MAX_COMPL_END_TH']     = float(self.primerMaxPairEndSpinBox.value())
            settings.params['PRIMER_MAX_HAIRPIN']               = float(self.primerMaxSelfHairpinSpinBox.value())
            settings.params['PRIMER_PAIR_MAX_HAIRPIN']          = float(self.primerMaxPairHairpinSpinBox.value())
            settings.params['PRIMER_MAX_END_STABILITY']         = float(self.primerMaxEndStabilitySpinBox.value())
            if self.primerGCClampCheckBox.isChecked:
                settings.params['PRIMER_GC_CLAMP']              = 1
            else:
                settings.params['PRIMER_GC_CLAMP']              = 0
            # create the primers table
            self.createPrimersTable()
                    
        if self.tagPrimersCheckBox.isChecked() or self.pigtailPrimersCheckBox.isChecked():
            # setup the settings for tagging primers
            tag_settings = primer.Settings()
            tag_settings.reduced(os.path.join(os.getcwd(), 'primer3_config/'), PRIMER_PICK_ANYWAY=1)
            # create the tagged primers table
            self.createTaggedPrimersTable()
            
        for record in SeqIO.parse(open(self.infile,'rU'), 'fasta'):
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            # add matches attribute to record
            record.matches = {}
            record.combined = {}
            record.primers = {}
            # search for the motif in the record.seq; store results in matches
            # attribute, biatch
            for motif in self.collection:
                self.searchForMotif(record, motif)
            if self.combineLociCheckBox.isChecked():
                record = self.combineLoci(record, int(self.combineLociDistanceSpinBox.value()))
            if self.designPrimersCheckBox.isChecked():
                if record.combined: 
                    matches = record.combined
                else:
                    matches = record.matches
                for match in matches:
                    primers = ()
                    for locations in matches[match]:
                        target = '%s,%s' % (locations[0][0], locations[0][1]-locations[0][0])
                        primer3 = primer.Primers(binary=os.path.join(os.getcwd(), 'primer3_core'))
                        primer3.pick(settings, sequence=str(record.seq), target=target, name = 'primers')
                        if primer3.primers_designed:
                            if self.pigtailPrimersCheckBox.isChecked() and not self.tagPrimersCheckBox.isChecked():
                                primer3.pigtail(tag_settings, str(self.pigtailPrimersTagLineEdit.text()))
                            elif self.tagPrimersCheckBox.isChecked():
                                cag, m13r, custom = None, None, None
                                if self.cagTagCheckBox.isChecked(): cag = 'CAGTCGGGCGTCATCA'
                                if self.m13rCheckBox.isChecked(): m13r = 'GGAAACAGCTATGACCAT'
                                if self.customTagCheckBox.isChecked(): custom = str(self.customTagLineEdit.text())
                                primer3.tag(tag_settings, CAG=cag, M13R=m13r, Custom = custom)
                                if primer3.tagged_good:
                                    if self.pigtailPrimersCheckBox.isChecked():
                                        primer3.pigtail(tag_settings, str(self.pigtailPrimersTagLineEdit.text()))
                        primers += ((primer3),)
                    record.primers[match] = primers
            # insert the referential integrity data first
            self.cur.execute('''INSERT INTO records (id, name) VALUES (?,?)'''\
                , (index, record.name))
            # if we are keeping the sequence reads pickle the sequence record, 
            # because we're gonna use it later and insert it into sqlite (we 
            # have to run Binary on it, first)
            if self.keepDbaseSequenceRecords.isChecked():
                rPickle = cPickle.dumps(record, 1)
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
            
            if self.combineLociCheckBox.isChecked():
                for match in record.combined:
                    for motif in record.combined[match]:
                        self.cur.execute('''INSERT INTO combined \
                            (records_id, id, motif, start, end, preceding, \
                            following, members) VALUES (?,?,?,?,?,?,?,?)''', \
                            (index, combine_index, match, motif[0][0], \
                            motif[0][1], motif[1], motif[2], motif[3]))
                        for m in motif[4]:
                            self.cur.execute('''INSERT INTO combined_components \
                            (records_id, combined_id, motif, length) VALUES \
                            (?,?,?,?)''', (index, combine_index, m[0], m[1]))
                        combine_index += 1
            
            if record.primers:
                if not self.combineLociCheckBox.isChecked():
                    self.storePrimers(index, msat_index, record)
                else:
                    self.storePrimers(index, combine_index, record)
            
            # this could be nested in if-then above
            if record.primers and (self.pigtailPrimersCheckBox.isChecked() or self.tagPrimersCheckBox.isChecked()):
                if not self.combineLociCheckBox.isChecked():
                    self.storeTaggedPrimers(index, msat_index, record)
                else:
                    self.storeTaggedPrimers(index, combine_index, record)
            
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
        # create the new primers table
        query = ('''CREATE TABLE primers (
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
            )''')
        self.cur.execute(query)
        self.conn.commit()
        
    def createTaggedPrimersTable(self):
        '''add tables to the dbase to hold the tagged primers'''
        # create the new primers table
        self.cur.execute('''CREATE TABLE tagged (
            records_id int,
            msats_id int,
            primer int,
            best int,
            tag text,
            tagged text,
            tag_seq test,
            common text,
            pigtail_tagged text,
            pigtail_tag_seq text,
            pigtail_common text,
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
    
    def storePrimers(self, record_id, msat_id, record):
        '''store primers in the database'''
        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        for motif, loci in record.primers.iteritems():
            for primers in loci:
                for i,p in primers.primers.iteritems():
                    if i != 'metadata':
                        # create a copy of the dict, to which we add the
                        # FOREIGN KEY reference
                        td = p.copy()
                        td['RECORD_ID'] = record_id
                        td['MSAT_ID']   = msat_id
                        td['PRIMER']    = i
                        query = ('''INSERT INTO primers VALUES (
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
                            )''')
                        self.cur.execute(query, td)
        self.conn.commit()
    
    def storeTaggedPrimers(self, record_id, msat_id, record, best=None):
        '''store primers in the database'''
        for motif, loci in record.primers.iteritems():
            for primers in loci:
                if primers.tagged_good:
                    for i,p in primers.tagged_good.iteritems():
                        if i != 'metadata':
                            # create a copy of the dict, to which we add the
                            # FOREIGN KEY reference
                            td = p.copy()
                            td['RECORD_ID'] = record_id
                            td['MSAT_ID']   = msat_id
                            td['BEST']      = 0
                            td['PRIMER'], td['TAG'], td['TAGGED'] = i.split('_')
                            self.cur.execute('''INSERT INTO tagged VALUES (
                                :RECORD_ID,
                                :MSAT_ID,
                                :PRIMER,
                                :BEST,
                                :TAG,
                                :PRIMER_TAGGED,
                                :PRIMER_TAG,
                                :PRIMER_TAG_COMMON_BASES,
                                :PRIMER_PIGTAILED,
                                :PRIMER_PIGTAIL_TAG,
                                :PRIMER_PIGTAIL_TAG_COMMON_BASES,
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
                    #QtCore.pyqtRemoveInputHook()
                    #pdb.set_trace()
                    if primers.tagged_best:
                        best = primers.tagged_best.keys()[0].split('_')
                        # we're using real side names now
                        if best[2]=='r':
                            side = 'RIGHT'
                        else:
                            side = 'LEFT'
                        self.cur.execute('''UPDATE tagged 
                            SET best = 1 WHERE 
                            records_id = ? 
                            AND msats_id = ?
                            AND primer = ? 
                            AND tag = ? 
                            AND tagged = ?''',
                            (record_id, msat_id, best[0], best[1], side))
        self.conn.commit()
        
    def outputWriter(self, out, extension, header, data):
        outpath = os.path.join(str(self.outdir), out)
        f = open(outpath, 'w')
        if extension == 'csv':
            h = ','.join(['\"{0}\"'.format(str(x)) for x in header])
            f.write('{0}\n'.format(h)) 
            for d in data:
                s = ','.join(['\"{0}\"'.format(str(x)) for x in d])
                f.write('{0}\n'.format(s))
        elif extension == 'tdt':
            h = '\t'.join([str(x) for x in header])
            f.write('{0}\n'.format(h))
            for d in data:
                s = '\t'.join([str(x) for x in d])
                f.write('{0}\n'.format(s))
        f.close()
    
    def outputResults(self):
        '''write results to some sort of output file'''
        # set the output format
        if self.commaSeparatedRadioButton.isChecked():
            extension = 'csv'
        elif self.tabDelimitedRadioButton.isChecked():
            extension = 'tdt'
        
        # if we want the actual, per-locus data returned
        if self.repeatsCheckBox.isChecked():
            out = 'msatcommander.microsatellites.%s' % extension
            self.cur.execute('''SELECT records.name, microsatellites.* from 
                records, microsatellites where 
                records.id = microsatellites.records_id''')
            header = [_[0] for _ in self.cur.description]
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
            
        # if we combined loci, and want those data output
        if self.combineLociCheckBox.isChecked():
            out = 'msatcommander.microsatellites.combined.%s' % extension
            self.cur.execute('''SELECT records.name, combined.* 
                from records, combined where 
                records.id = combined.records_id''')
            header = [_[0] for _ in self.cur.description]
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
        
        # if we designed primers, and we didn't pigtail them, and we want the
        # untagged primer data
        #
        # This will return the best primer designed for each locus 
        # (combined or not)
        if self.primersCheckBox.isChecked() \
                and not self.tagPrimersCheckBox.isChecked() \
                and not self.pigtailPrimersCheckBox.isChecked():
            out = 'msatcommander.primers.%s' % extension
            # get the best unlabelled primers
            self.cur.execute('''SELECT records.name, primers.* from 
                records, primers where records.id = primers.records_id and 
                primers.primer = 0''')
            header = [_[0] for _ in self.cur.description]
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
        
        # if we designed primers, and we pigtailed primers, and we want the
        # untagged primer data
        #
        # This will return *only* the primers matching the best pigtailed
        # primers
        if self.designPrimersCheckBox.isChecked() \
                and self.primersCheckBox.isChecked() \
                and self.pigtailPrimersCheckBox.isChecked() \
                and not self.tagPrimersCheckBox.isChecked():
            out = 'msatcommander.primers.%s' % extension
            # get best primers based on best pigtailed primers
            query = '''SELECT 
                records.name,
                primers.records_id,
                primers.msats_id,
                primers.primer,
                primers.left,
                tagged.tag,
                tagged.tagged,
                tagged.tag_seq,
                tagged.left_sequence,
                primers.left_tm,
                primers.left_gc,
                primers.left_self_end,
                primers.left_self_any,
                primers.left_hairpin,
                primers.left_end_stability,
                primers.left_penalty,
                primers.right,
                tagged.right_sequence,
                primers.right_tm,
                primers.right_gc,
                primers.right_self_end,
                primers.right_self_any,
                primers.right_hairpin,
                primers.right_end_stability,
                primers.right_penalty,
                primers.pair_product_size,
                primers.pair_compl_end,
                primers.pair_compl_any,
                primers.pair_penalty
                FROM records, primers, tagged WHERE 
                primers.records_id = records.id 
                AND primers.records_id = tagged.records_id 
                AND primers.msats_id = tagged.msats_id
                AND primers.primer = tagged.primer
                AND tagged.best = 1'''
            self.cur.execute(query)
            #QtCore.pyqtRemoveInputHook()
            #pdb.set_trace()
            header = [_[0] for _ in self.cur.description]
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
        
        # if we designed primers and we want to ouput tagged primers - they
        # may or may not be pigtailed
        #
        # This will return the best tagged primers and only the 
        # untagged primers matching the best tagged primers
        if self.designPrimersCheckBox.isChecked() \
                and self.taggedPrimersCheckBox.isChecked():
            out = 'msatcommander.primers.%s' % extension
            # get best primers based on best tagged primers
            query = '''SELECT records.name, primers.* FROM records, primers, 
                tagged WHERE 
                primers.records_id = records.id
                AND primers.records_id = tagged.records_id
                AND primers.msats_id = tagged.msats_id
                AND primers.primer = tagged.primer
                AND tagged.best = 1'''
            self.cur.execute(query)
            header = [_[0] for _ in self.cur.description]
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)
            # get the tagged primers
            out = 'msatcommander.tagged.%s' % extension
            self.cur.execute('''SELECT records.name, tagged.* FROM records,
                tagged WHERE 
                tagged.records_id = records.id 
                AND tagged.best = 1''')
            header = [_[0] for _ in self.cur.description]
            data = self.cur.fetchall()
            self.outputWriter(out, extension, header, data)

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())