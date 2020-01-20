#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time, sys
import argparse, sys, os
from PyQt5.QtWidgets import * # (QWidget, QLabel, QLineEdit, QTextEdit, QGridLayout, QApplication, QPushButton)


from PyQt5.QtGui import QColor

sys.setrecursionlimit(2010)

whitecolor = QColor(255,255,255)
blackColor = QColor(0, 0, 0)
blueColor = QColor(0, 0, 255)
greenColor = QColor(0, 150, 0)
redColor = QColor(255, 0, 0)

class TermAnnotator(QWidget):
    
	def __init__(self, sfile, tfile, termstfile, termtsfile): #, nosyn_check):
		super().__init__()
	
		self.bicor = {}
		self.sent_ID = 1

		self.termTableST = None
		self.max_ann_sterm_len_ST = 0 
		self.termTableTS = None
		self.max_ann_sterm_len_TS = 0

		self.mono_syn_S = {}
		self.mono_syn_T = {}

		self.mono_syn_S_depth_2 = {}
		self.mono_syn_T_depth_2 = {}

		self.TermFinderTableStoT = None 
		self.max_sterm_len_ST = 0
		self.TermFinderTableTtoS = None
		self.max_sterm_len_TS = 0

		self.termsugg_start_s = []
		self.termsugg_end_s = []
		self.termsugg_start_t = []
		self.termsugg_end_t = []
		self.currentsugg = 0


		self.termsugg_mono_S = {}
		self.termsugg_mono_T = {}
		self.currentsugg_mS = 0
		self.currentsugg_mT = 0
		self.currentsynlist = []
		self.currentsyn_m_pos = 0
		self.m_Flag = 0

		self.accepted_mS = {}
		self.accepted_mT = {}
		# current sentence ->
		self.curr_rejected_mS = {}
		self.curr_rejected_mT = {}

		self.accepted_start_s = []
		self.accepted_end_s = []
		self.accepted_start_t = []
		self.accepted_end_t = []

		self.rejected_start_s = []
		self.rejected_end_s = []
		self.rejected_start_t = []
		self.rejected_end_t = []

		self.rejected_pairs = {}
		self.rejected_mono_s = {}
		self.rejected_mono_t = {}
		
		#self.annotated_sent_pairs = {}
		self.annotated_sent_ID = set()

		self.starttime = time.time()
		self.monostarttime = time.time()
		self.endtime = time.time()
		self.no_acc_tp = 0
		self.no_rej_tp = 0
		self.no_skip_tp = 0
		self.no_new_tp = 0

		self.no_acc_m1 = 0
		self.no_rej_m1 = 0

		self.no_acc_m2 = 0
		self.no_rej_m2 = 0
		self.no_skip_m = 0

		self.syndic_s ={} 
		self.syndic_t = {}

		self.readfiles(sfile, tfile)
		if termstfile:
			(self.TermFinderTableStoT, self.max_sterm_len_ST) = self.read_term_file(termstfile)
		if termtsfile:
			(self.TermFinderTableTtoS, self.max_sterm_len_TS) = self.read_term_file(termtsfile)
		
		# Uses already annotated term-pairs to annotate new sentence-pairs ### also uses rejected terms --> this speed up the process
		self.load_annotated_terms()
		# Uses already annotated mono term with its synonym suggestions # also uses rejected suggs --> this also speed up 
		self.load_monolingual_synonymfiles()
		#loads ids of sent-pairs that have been already annotaed / rejected for spelling error or other errors
		self.load_annotated_sent_pairs()

		# For colour term if from already annotated
		self.from_annot = []

		#self.nosyn_check = nosyn_check
		self.initUI()

	def load_monolingual_synonymfiles(self):
		if not os.path.isfile("syn.src") or not os.path.isfile("syn.trg"):
			return 0

		fsyns = open("syn.src")
		fsynt = open("syn.trg")

		for line in fsyns:
			synlist = line.strip().split('\t')
			s1 = synlist[0].lower()
			s2 = synlist[1].lower()

			if s1 in self.mono_syn_S:
				if s2 not in self.mono_syn_S[s1]:
					self.mono_syn_S[s1].append(s2)
			else:
				self.mono_syn_S[s1] = [s2]
			if s2 in self.mono_syn_S:
				if s1 not in self.mono_syn_S[s2]:
					self.mono_syn_S[s2].append(s1)
			else:
				self.mono_syn_S[s2] = [s1]

		for line in fsynt:
			synlist = line.strip().split('\t')
			#print (synlist)
			t1 = synlist[0].lower()
			t2 = synlist[1].lower()

			if t1 in self.mono_syn_T:
				if t2 not in self.mono_syn_T[t1]:
					self.mono_syn_T[t1].append(t2)
			else:
				self.mono_syn_T[t1] = [t2]
			if t2 in self.mono_syn_T:
				if t1 not in self.mono_syn_T[t2]:
					self.mono_syn_T[t2].append(t1)
			else:
				self.mono_syn_T[t2] = [t1]

		fsyns.close()
		fsynt.close()

		if not os.path.isfile("syn_rejected.src") or not os.path.isfile("syn_rejected.trg"):
			return 0

		fsyns_rejected = open("syn_rejected.src")
		fsynt_rejected = open("syn_rejected.trg")

		for line in fsyns_rejected:
			linelist = line.lower().strip().split('\t')
			self.rejected_mono_s[(linelist[0], linelist[1])] = 1
		
		for line in fsynt_rejected:
			linelist = line.lower().strip().split('\t')
			self.rejected_mono_t[(linelist[0], linelist[1])] = 1
			
		fsyns_rejected.close()
		fsynt_rejected.close()

	def cleanterm(self, term):

		term = " ".join(term.split())
		return term.strip()

	def duplicate_eliminate(self, filename):
		if not os.path.isfile(filename):
			return 0

		lineset = set()
		f = open(filename)
		for line in f:
			line = line.lower().strip()
			if '\t' in line:
				linel = line.split('\t')
				tline = ""
				for term in linel:
					tline+=self.cleanterm(term)+"\t"
				lineset.add(tline.strip())
			#else:
				#lineset.add(line) #
		f.close()

		fw = open(filename, "w")
		for line in sorted(lineset):
			fw.write(line+"\n")
		fw.close()
		
	def remove_from_rejected(self, correctfile, rejectedfile):
		if not os.path.isfile(correctfile) or not os.path.isfile(rejectedfile):
			return 0

		lineset = set()
		f = open(correctfile)
		for line in f:
			lineset.add(line)
		f.close()

		lineset2 = set()
		f = open(rejectedfile)
		for line in f:
			lineset2.add(line)
		f.close()

		fw = open(rejectedfile, "w")
		for line in sorted(lineset2):
			if line not in lineset:
				fw.write(line)
		fw.close()
		
		
	

	def eliminate_duplicates_from_term_bases(self):

		self.duplicate_eliminate("syn.src")
		self.duplicate_eliminate("syn.trg")
		self.duplicate_eliminate("syn_rejected.src")
		self.duplicate_eliminate("syn_rejected.trg")
		self.duplicate_eliminate("annotated_terms")
		self.duplicate_eliminate("rejected_term_pairs")

		self.remove_from_rejected("annotated_terms","rejected_term_pairs")
		self.remove_from_rejected("syn.src","syn_rejected.src")
		self.remove_from_rejected("syn.trg","syn_rejected.trg")

	def load_annotated_sent_pairs(self):
		#if not os.path.isfile("anncorpus.src") or not os.path.isfile("anncorpus.trg") or not os.path.isfile("corpus.src") or not os.path.isfile("corpus.trg") or not os.path.isfile("annID.txt"):
		if os.path.isfile("annID.txt"):
			fid = open("annID.txt")
			for idstr in fid:
				self.annotated_sent_ID.add(int(idstr.strip()))
			fid.close()

		if os.path.isfile("omittedID.txt"):
			fid = open("omittedID.txt")
			for idstr in fid:
				self.annotated_sent_ID.add(int(idstr.strip()))
			fid.close()

			

		#fs = open("anncorpus.src")
		#ft = open("anncorpus.trg")
		#fsOrig = open("corpus.src")
		#ftOrig = open("corpus.trg")

		#for tagged_s_line in fs:
		#	tagged_t_line = ft.readline()
		#	s_line = fsOrig.readline()
		#	t_line = ftOrig.readline()
		#	#self.annotated_sent_pairs[(s_line.lower(), t_line.lower())] = (tagged_s_line, tagged_t_line)

		#fid.close()
		#ft.close()
		#fs.close()
		#fsOrig.close()
		#ftOrig.close()
	
	def readfiles(self, sfile, tfile):
		fs = open(sfile)
		ft = open(tfile)
		c =  0
		for lines in fs:
			lines = lines.strip()
			linet = ft.readline().strip()
			self.bicor[c] = (lines, linet) 
			c+=1
		fs.close()
		ft.close() 
	
	def read_term_file(self, termfile): 
		tf = open(termfile)
		table = []
		max_sterm_len = 0
		for term in tf:
			termlist = term.lower().strip().split("\t")
			sterm = termlist[0]
			length_s = len(sterm)
			if len(sterm.split()) > max_sterm_len:
				max_sterm_len = len(sterm.split())
			while len(table) < length_s+1:
				table.append({})
			if not sterm in table[length_s]:
				table[length_s][sterm] = termlist[1:]
			else:
				for tterm in termlist[1:]:
					if tterm not in table[length_s][sterm]:
						table[length_s][sterm].append(tterm) #extend(termlist[1:])
		tf.close()

		return (table, max_sterm_len)

	def load_annotated_terms(self):
		if not os.path.isfile("annotated_terms"):
			return 0
	
		max_sterm_len = 0
		max_tterm_len = 0	

		termTableST = []
		termTableTS = []

		fw = open("annotated_terms")
		for line in fw:
			st = line.lower().strip().split("\t")
			sterm = st[0]
			tterm = st[1]
			length_s = len(sterm)
			length_t = len(tterm)
			if len(sterm.split()) > max_sterm_len:
				max_sterm_len = len(sterm.split())
			if len(tterm.split()) > max_tterm_len:
				max_tterm_len = len(tterm.split())

			while len(termTableST) < length_s+1:
				termTableST.append({})
			while len(termTableTS) < length_t+1:
				termTableTS.append({})

			if not sterm in termTableST[length_s]:
				termTableST[length_s][sterm] = [tterm]
			else:
				termTableST[length_s][sterm].append(tterm)

			if not tterm in termTableTS[length_t]:
				termTableTS[length_t][tterm] = [sterm]
			else:
				termTableTS[length_t][tterm].append(sterm)
			
		fw.close()
		self.termTableST = termTableST
		self.max_ann_sterm_len_ST = max_sterm_len
		self.termTableTS = termTableTS
		self.max_ann_sterm_len_TS = max_tterm_len

		if not os.path.isfile("rejected_term_pairs"):
			return 0

		fw = open("rejected_term_pairs")
		for line in fw:
			st = line.strip().lower().split("\t")
			sterm = st[0]
			tterm = st[1]

			self.rejected_pairs[(sterm,tterm)] = 1

		fw.close()
				

	def save_this_annotated_sent_pairs(self):
		self.endtime = time.time()
		self.save_sent_pairs(self.ssent_O_Edit.toPlainText(), self.tsent_O_Edit.toPlainText(), self.bicor[self.sent_ID-1][0], self.bicor[self.sent_ID-1][1], self.ssent_al_Edit.toPlainText(), self.tsent_al_Edit.toPlainText(), self.sent_ID)
		self.update_term_table()
		#print (self.endtime - self.starttime)
		# Load available sentence starting from self.sent_ID
		self.sent_ID+=1
		self.get_next_available_ID()
		self.loadSentence()

	def save_sent_pairs(self, ann_syn_src, ann_syn_trg, sline, tline, al_s, al_t, sent_ID):
		fws = open("anncorpus.src", "a")
		fwt = open("anncorpus.trg", "a")
		#fwsOrig = open("corpus.src", "a")
		#fwtOrig = open("corpus.trg", "a")
		fwid = open("annID.txt", "a")
		timef = open("timestat.txt", "a")

		ann_syn_src = ann_syn_src.replace("\n", "|")
		ann_syn_trg = ann_syn_trg.replace("\n", "|")
	
		fws.write(sline + " ||| "+ann_syn_src.replace("\n", "|")+" ||| "+al_s+"\n")
		fwt.write(tline + " ||| "+ann_syn_trg.replace("\n", "|")+" ||| "+al_t+"\n")
		#fwsOrig.write(sline+"\n")
		#fwtOrig.write(tline+"\n")

		totmontime = str(self.monostarttime - self.starttime)
		if len(self.termsugg_mono_S) ==0:
			totmontime = "0.0"
		tottime = str(self.endtime - self.starttime)
		fwid.write(str(sent_ID)+"\n")
		#print (str(self.monostarttime), str(self.starttime), str(self.endtime))
		#print (totmontime, tottime)
		#timef.write(str(sent_ID)+"\t"+str(self.monostarttime - self.starttime)+"\t"+str(self.endtime - self.starttime)+"\t")
		timef.write(str(sent_ID)+"\t"+totmontime+"\t"+tottime+"\t")
		timef.write(str(self.no_acc_tp)+"\t"+str(self.no_rej_tp)+"\t"+str(self.no_skip_tp)+"\t"+str(self.no_new_tp)+"\t"+str(self.no_acc_m1)+"\t"+str(self.no_rej_m1)+"\t"+str(self.no_acc_m2)+"\t"+str(self.no_rej_m2)+"\t"+str(self.no_skip_m)+"\n")

		timef.close()
		fwid.close()
		fwt.close()
		fws.close()
		#fwsOrig.close()
		#fwtOrig.close()

		self.write_deleted_syn("deleted_syn.src", "syn.src", al_s, ann_syn_src, self.syndic_s, self.sent_ID)
		self.write_deleted_syn("deleted_syn.trg", "syn.trg", al_t, ann_syn_trg, self.syndic_t, self.sent_ID)

		#self.annotated_sent_pairs[(sline.lower(), tline.lower())] = (ann_src, ann_trg)
		self.annotated_sent_ID.add(sent_ID)

	def return_set(self, synonym):
		s = set()
		for syn in synonym.split(":"):
			syn = syn.strip()
			if syn:
				s.add(syn)
		return s

	def addsyn(self,term, syn, mono_syn):
		
		if term.lower() in mono_syn:
			if syn.lower() not in mono_syn[term.lower()]:
				mono_syn[term.lower()].append(syn.lower())
		else:
			mono_syn[term.lower()] = [syn.lower()]

		if syn.lower() in mono_syn:
			if term.lower() not in mono_syn[syn.lower()]:
				mono_syn[syn.lower()].append(term.lower())
		else:
			mono_syn[syn.lower()] = [term.lower()]

	def write_deleted_syn(self, filename, synfile, alignment, curr_syn, syndic, sentid):

		fwdel = open(filename, "a")
		fwsyn = open(synfile, "a")


		aldic = {}
		for al in alignment.strip().split():
			aldic[al.split(":")[0]] = al.split(":")[1]

		monosyndictionary = self.mono_syn_S
		if synfile == "syn.trg":
			monosyndictionary = self.mono_syn_T

		for eachsyn in curr_syn.split("|"):
			synlist = eachsyn.split(":")
			edited_syn_set = self.return_set(':'.join(synlist[2:]).strip(":"))

			print ("edited set", edited_syn_set)
			for syn in edited_syn_set:
				term = synlist[1].strip().strip(":")
				fwsyn.write(term+"\t"+syn+"\n")
				self.addsyn(term, syn, monosyndictionary)

			if len(synlist) >1:
				key = synlist[0]+":"+synlist[1]+":"
				if key in syndic:
					old_syn_set = syndic[key]
					print ("old synset", old_syn_set)
					for syn in old_syn_set:
						if syn not in edited_syn_set:
							key = synlist[0].replace("trm","")
							fwdel.write(str(sentid)+"\t"+synlist[0]+"\t"+aldic[key]+"\t"+syn+"\n")
		fwdel.close()
		fwsyn.close()



	def update_term_table(self):

		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		fsyns = open("syn.src", "a")
		fsynt = open("syn.trg", "a")
		fw = open("annotated_terms", "a")
		fa_syns = open("syn_rejected.src", "a")
		fa_synt = open("syn_rejected.trg", "a")


		for i in range (len(self.accepted_start_s)):
			sterm = ' '.join(sline.split()[self.accepted_start_s[i]:self.accepted_end_s[i]+1]).strip()
			tterm = ' '.join(tline.split()[self.accepted_start_t[i]:self.accepted_end_t[i]+1]).strip()
		
			sterm = sterm.lower()
			tterm = tterm.lower()			

			print ("writing", sterm, tterm)	
			fw.write(sterm+"\t"+tterm+"\n")

			if len(sterm.split()) > self.max_ann_sterm_len_ST:
				self.max_ann_sterm_len_ST = len(sterm.split())
			length_s = len(sterm)
			if self.termTableST:
				while len(self.termTableST) < length_s+1:
					self.termTableST.append({})
			if not sterm in self.termTableST[length_s]:
				self.termTableST[length_s][sterm] = [tterm]
			else:
				if tterm not in self.termTableST[length_s][sterm]:
					self.termTableST[length_s][sterm].append(tterm)

			if len(tterm.split()) > self.max_ann_sterm_len_TS:
				self.max_ann_sterm_len_TS = len(tterm.split())

			length_t = len(tterm)
			if self.termTableTS:
				while len(self.termTableTS) < length_t+1:
					self.termTableTS.append({})
			if not tterm in self.termTableTS[length_t]:
				self.termTableTS[length_t][tterm] = [sterm]
			else:
				if sterm not in self.termTableTS[length_t][tterm]:
					self.termTableTS[length_t][tterm].append(sterm)

			###########################################################
			# Saving Monolingual synonyms for sterm and tterm
			
			key_s = (self.accepted_start_s[i], self.accepted_end_s[i])
			#if key_s in self.accepted_mS:
			#	accepted_syn = self.accepted_mS[key_s]
			#	for syn in accepted_syn:
					#fsyns.write(sterm.lower()+"\t"+syn.lower()+"\n")
			#		if sterm.lower() in self.mono_syn_S:
			#			if syn.lower() not in self.mono_syn_S[sterm.lower()]:
			#				self.mono_syn_S[sterm.lower()].append(syn.lower())
			#		else:
			#			self.mono_syn_S[sterm.lower()] = [syn.lower()]
			#		if syn.lower() in self.mono_syn_S:
			#			if sterm.lower() not in self.mono_syn_S[syn.lower()]:
			#				self.mono_syn_S[syn.lower()].append(sterm.lower())
			#		else:
			#			self.mono_syn_S[syn.lower()] = [sterm.lower()]
					
			key_t = (self.accepted_start_t[i], self.accepted_end_t[i])
			#if key_t in self.accepted_mT:
			#	accepted_syn = self.accepted_mT[key_t]
			#	for syn in accepted_syn:
					#fsynt.write(tterm.lower()+"\t"+syn.lower()+"\n")
			#		if tterm.lower() in self.mono_syn_T:
			#			if syn.lower() not in self.mono_syn_T[tterm.lower()]:
			#				self.mono_syn_T[tterm.lower()].append(syn.lower())
			#		else:
			#			self.mono_syn_T[tterm.lower()] = [syn.lower()]
			#		if syn.lower() in self.mono_syn_T:
			#			if tterm.lower() not in self.mono_syn_T[syn.lower()]:
			#				self.mono_syn_T[syn.lower()].append(tterm.lower())
			#		else:
			#			self.mono_syn_T[syn.lower()] = [tterm.lower()]

			if key_s in self.curr_rejected_mS:
				rejected_syn = self.curr_rejected_mS [key_s]
				for rsyn in rejected_syn:
					fa_syns.write(sterm+"\t"+rsyn.lower()+"\n")
					self.rejected_mono_s[(sterm, rsyn.lower())] = 1
			if key_t in self.curr_rejected_mT:
				rejected_syn = self.curr_rejected_mT [key_t]
				for rsyn in rejected_syn:
					fa_synt.write(tterm+"\t"+rsyn.lower()+"\n")
					self.rejected_mono_t[(tterm, rsyn.lower())] = 1
					
				

		fwreject = open("rejected_term_pairs", "a")
		for i in range (len(self.rejected_start_s)):
			sterm = ' '.join(sline.split()[self.rejected_start_s[i]:self.rejected_end_s[i]+1]).strip().lower()
			tterm = ' '.join(tline.split()[self.rejected_start_t[i]:self.rejected_end_t[i]+1]).strip().lower()

			fwreject.write(sterm+"\t"+tterm+"\n")
			self.rejected_pairs[(sterm, tterm)] = 1

		fa_skip = open("skipped_term_pairs", "a")
		for i in range(len(self.termsugg_start_s)):
				s_s = self.termsugg_start_s[i]
				e_s = self.termsugg_end_s[i]
				s_t = self.termsugg_start_t[i]
				e_t = self.termsugg_end_t[i]

				accepted = False
				for x in range (len(self.accepted_start_s)):
					if s_s == self.accepted_start_s[x] and e_s == self.accepted_end_s[x] and s_t == self.accepted_start_t[x] and e_t == self.accepted_end_t[x]:
						accepted = True
						break

				if not accepted:
					sterm = ' '.join(sline.split()[s_s: e_s+1]).strip().lower()
					tterm = ' '.join(tline.split()[s_t: e_t+1]).strip().lower()
					fa_skip.write(str(self.sent_ID)+"\t"+ str(s_s)+"\t"+ str(e_s)+"\t"+ str(s_t)+"\t"+ str(e_t)+"\t"+sterm+"\t"+tterm+"\n")

		

		#fa_syns = open("syn_rejected.src", "a")
		#fa_synt = open("syn_rejected.trg", "a")

	
		#for (terms, syns) in self.curr_rejected_mS:
		#	self.rejected_mono_s[(terms.lower(), syns.lower())] = 1
		#	fa_syns.write(terms.lower()+"\t"+syns.lower()+"\n")
		#for (termt, synt) in self.curr_rejected_mT:
		#	self.rejected_mono_t[(termt.lower(), synt.lower())] = 1
		#	fa_synt.write(termt.lower()+"\t"+synt.lower()+"\n")
		
		fa_syns.close()
		fa_synt.close()		
		fwreject.close()
		fw.close()
		fsyns.close()
		fsynt.close()		
		fa_skip.close()

	def displayNextSentence(self):
		self.sent_ID+=1
		self.get_next_available_ID()
		self.loadSentence()

	def omit_this_sent(self, sent_ID):
		fwid = open("omittedID.txt", "a")
		fwid.write(str(sent_ID)+"\n")
		fwid.close()	

	def RejectThisSent_displayNextSentence(self):
		self.omit_this_sent(self.sent_ID)
		self.sent_ID+=1
		self.get_next_available_ID()
		self.loadSentence()
		

	def getSentenceID(self):
		self.sent_ID = int(self.textbox.text())
		self.get_next_available_ID()
		self.loadSentence()

	def get_next_available_ID(self):
		if len(self.annotated_sent_ID) == len(self.bicor):
			print ("Term Annotation Task Done!")
			sys.exit(0)
		if self.check_if_already_annotated():
			self.sent_ID+=1
			self.get_next_available_ID()
			if self.sent_ID > len(self.bicor):
				# Search again from first sentence
				self.sent_ID = 1

	def loadSentence(self):
		self.textbox.setText(str(self.sent_ID))
		print ("Sentence ID:"+ str(self.sent_ID))
		self.initialise_members()	
		self.get_suggested_annoatations()
		self.get_additional_sugg_from_TermFinder()
		self.remove_rejected_suggestions()
		self.starttime = time.time()
		self.monostarttime = time.time()	
		self.start_annotation()
		

	def remove_rejected_suggestions(self):
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		s_s = []
		s_e = []
		t_s = []
		t_e = []
		for i in range(len(self.termsugg_start_s)):
			sterm = ' '.join(sline.split()[self.termsugg_start_s[i]: self.termsugg_end_s[i]+1]).strip().lower()
			tterm = ' '.join(tline.split()[self.termsugg_start_t[i]: self.termsugg_end_t[i]+1]).strip().lower()
			#print (sterm, tterm)
			if (sterm, tterm) not in self.rejected_pairs:
				#print (sterm, tterm)
				s_s.append(self.termsugg_start_s[i])
				s_e.append(self.termsugg_end_s[i])
				t_s.append(self.termsugg_start_t[i])
				t_e.append(self.termsugg_end_t[i])

		self.termsugg_start_s = s_s
		self.termsugg_end_s = s_e
		self.termsugg_start_t = t_s
		self.termsugg_end_t = t_e

	def check_if_already_annotated(self):
		#print (self.sent_ID)
		#print (self.annotated_sent_ID)
		if self.sent_ID in self.annotated_sent_ID:
			return True
		#if (self.bicor[self.sent_ID-1][0].lower(), self.bicor[self.sent_ID-1][1]) in self.annotated_sent_pairs:
		#	tup_ann_pair = self.annotated_sent_pairs[(self.bicor[self.sent_ID-1][0].lower(), self.bicor[self.sent_ID-1][1])]
			# Update ann_cor, cor, and ID files; just updated ann_sent_pairs dic and sent_id set
		#	self.save_sent_pairs(tup_ann_pair[0], tup_ann_pair[1], self.bicor[sent_ID-1][0], self.bicor[self.sent_ID-1][1], self.sent_ID)
		#	return True
		return False

	def initialise_members(self):
		
		self.termsugg_mono_S = {}
		self.termsugg_mono_T = {}

		self.currentsugg_mS = 0
		self.currentsugg_mT = 0
		self.currentsyn_m_pos = 0
		self.currentsynlist = []
		self.accepted_mT = {}	
		self.accepted_mS = {}
		self.m_Flag = 0
		self.curr_rejected_mS = {}
		self.curr_rejected_mT = {}

		self.syndic_s  = {}
		self.syndic_t = {}

		######################################
		self.termsugg_start_s = []
		self.termsugg_end_s = []
		self.termsugg_start_t = []
		self.termsugg_end_t = []
		self.currentsugg = 0

		self.accepted_start_s = []
		self.accepted_end_s = []
		self.accepted_start_t = []
		self.accepted_end_t = []

		self.rejected_start_s = []
		self.rejected_end_s = []
		self.rejected_start_t = []
		self.rejected_end_t = []
	
		self.no_acc_tp = 0
		self.no_rej_tp = 0
		self.no_skip_tp = 0
		self.no_new_tp = 0

		self.from_annot = []

		#if not self.nosyn_check:	
		self.textbox_term.setText("")
		self.textbox_syn.setText("")

		#initialise monolingual members
		self.init_params()


	def start_annotation(self):
		retval = self.set_text()
		#print ("Return value from set_text: " + str(retval))
		#print ("Current Term index: "+str(self.currentsugg))
		self.displaysugg(retval)
		self.displayOutput()


	def set_text(self):
		totsugg = len(self.termsugg_start_s)
		cind = self.currentsugg
		if totsugg > 0 and cind < totsugg:
			text = str(cind+1) +"/"+str(totsugg)
			self.textbox_sugg1.setText(text)
			return 1
		elif cind >= totsugg:
			self.textbox_sugg1.setText("END")
			return 0
		else:
			self.textbox_sugg1.setText("NIL")
			return 0

	def displaysugg(self, retval):
		
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]
		cind = self.currentsugg
	
		if retval == 1:
			t= (self.termsugg_start_s [cind], self.termsugg_end_s [cind], self.termsugg_start_t [cind], self.termsugg_end_t [cind])
			colour = 0
			if t in self.from_annot: 
				self.textbox_sugg1.setStyleSheet("color: green;") #rgb(255, 0, 255);")
			else:
				colour = 1
				self.textbox_sugg1.setStyleSheet("color: blue;") #rgb(255, 0, 255);")

			#sugg_source = 
			self.get_line_with_term(sline, self.ssent_Edit, self.termsugg_start_s [cind], self.termsugg_end_s [cind], cind, colour)
			#sugg_target = 
			self.get_line_with_term(tline, self.tsent_Edit, self.termsugg_start_t [cind], self.termsugg_end_t [cind], cind, colour)
			#self.ssent_Edit.setText (sugg_source)
			#self.tsent_Edit.setText (sugg_target)
		else:
			self.ssent_Edit.setTextColor(blackColor)
			self.tsent_Edit.setTextColor(blackColor)
			self.ssent_Edit.setText (sline)
			self.tsent_Edit.setText (tline)
			
			self.textbox_sugg1.setStyleSheet("color: black;")

	def get_line_with_term(self, line, editbox, indexStart, indexEnd, termID, colour):
		global blackColor, greenColor, blueColor
		
		editbox.setText ("")
		editbox.setTextColor(blackColor)

		for index, word in enumerate (line.split()):
			if index == indexStart:
				#editbox.setBold(True)
				editbox.insertPlainText("<trm"+str (termID+1)+"> ")
				if colour == 0:
					editbox.setTextColor(greenColor)
				if colour == 1:
					editbox.setTextColor(blueColor)
			editbox.insertPlainText(word+ " ")		
			if index == indexEnd:
				#editbox.setBold(False)
				editbox.setTextColor(blackColor)
				editbox.insertPlainText("</trm> ")

	def retstr(self, lst):
		str1=""
		for syn in lst:
			syn = syn.strip()
			if syn != '-':
				str1+=syn+":"

		return str1.strip().strip(":")	


	def print_line(self, sent_O_Edit, line, accepted_start, accepted_end):
		global blackColor, greenColor, blueColor, redColor, whitecolor

		sent_O_Edit.setText("")
		for index, word in enumerate(line.split()):
			if index in accepted_start:
				indices = []
				Tlen = {}
				for i, x in enumerate(accepted_start):
					if x == index:
						term_id = i
						term_len = accepted_end[i] - accepted_start[i]
						Tlen[term_len] = term_id

				for key in sorted(Tlen.keys(), reverse=True):
                                        indices.append(Tlen[key])
				for termID in indices: #[::-1]:
					sent_O_Edit.setTextColor(redColor)
                                        #sent_O_Edit.insertPlainText("<trm"+str(termID+1)+"> ")

			sent_O_Edit.insertPlainText(word+ " ")
			if index in accepted_end:
				indices = [i for i, x in enumerate(accepted_end) if x == index]
				sent_O_Edit.setTextColor(blackColor)
                                #sent_O_Edit.setTextBackgroundColor(whitecolor)

                                #for termID in indices:




	def start_editing_syn_sugg(self):

		self.displayOutput_2()

	def displayOutput_2(self):
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		self.ssent_O_Edit.setText("")
		self.tsent_O_Edit.setText("")
		self.ssent_al_Edit.setText("")
		self.tsent_al_Edit.setText("")
		#self.ssent_O_Edit.insertPlainText(sline)
		#self.tsent_O_Edit.insertPlainText(tline)
		self.print_line(self.ssent_Edit, sline, self.accepted_start_s, self.accepted_end_s)
		self.print_line(self.tsent_Edit, tline, self.accepted_start_t, self.accepted_end_t)

		self.syndic_s  = {}
		self.syndic_t = {}


		IDs = ""
		IDt = ""		
		if len(self.accepted_start_s) > 0:
			for i in range (len(self.accepted_start_s)):
				sterm = ' '.join(sline.split()[self.accepted_start_s[i]:self.accepted_end_s[i]+1]).strip().lower()
				tterm = ' '.join(tline.split()[self.accepted_start_t[i]:self.accepted_end_t[i]+1]).strip().lower()
				#if sterm in self.mono_syn_S:
				key = (self.accepted_start_s[i], self.accepted_end_s[i])
				
				if key in self.accepted_mS: #self.termsugg_mono_S:
					self.ssent_O_Edit.insertPlainText("trm"+str(i+1)+":"+sterm+":"+ self.retstr(self.accepted_mS[key])+":\n")
					strkeyS = "trm"+str(i+1)+":"+sterm+":"
					self.syndic_s [strkeyS] = self.return_set(self.retstr(self.accepted_mS[key]))
				else:
					self.ssent_O_Edit.insertPlainText("trm"+str(i+1)+":"+sterm+":\n")

				key = (self.accepted_start_t[i], self.accepted_end_t[i])
				if key in self.accepted_mT:
					self.tsent_O_Edit.insertPlainText("trm"+str(i+1)+":"+tterm+":"+self.retstr(self.accepted_mT[key])+":\n")
					strkeyT = "trm"+str(i+1)+":"+tterm+":"
					self.syndic_t [strkeyT] = self.return_set(self.retstr(self.accepted_mT[key]))

				else:
					self.tsent_O_Edit.insertPlainText("trm"+str(i+1)+":"+tterm+":\n")

				#Alignment
				IDs += str(i+1)+":"+str(self.accepted_start_s[i]) +"-"+ str(self.accepted_end_s[i]) +" "
				IDt += str(i+1)+":"+str(self.accepted_start_t[i]) +"-"+ str(self.accepted_end_t[i]) +" "

			self.ssent_al_Edit.insertPlainText(IDs)
			self.tsent_al_Edit.insertPlainText(IDt)

		#else:
		#	self.ssent_O_Edit.insertPlainText(" ||| |||")
		#	self.tsent_O_Edit.insertPlainText(" ||| |||")


	def displayOutput(self):
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		if len(self.accepted_start_s) > 0:
			#out_source = self.get_line_with_accepted_term_pairs(  sline, self.accepted_start_s, self.accepted_end_s,self.accepted_mS, self.m_Flag, if self.currentsugg_mS)
			#out_target = self.get_line_with_accepted_term_pairs(tline, self.accepted_start_t, self.accepted_end_t,self.accepted_mT, self.m_Flag, if self.currentsugg_mT)
			currentsugg_mS = self.currentsugg_mS
			currentsugg_mT = self.currentsugg_mT
			if self.m_Flag == 2:
				currentsugg_mS = -1
			if self.m_Flag == 1:
                                currentsugg_mT = -1

			self.get_line_with_accepted_term_pairs (self.ssent_O_Edit, self.ssent_al_Edit, sline, self.accepted_start_s, self.accepted_end_s,self.accepted_mS, self.m_Flag, currentsugg_mS)
			self.get_line_with_accepted_term_pairs (self.tsent_O_Edit, self.tsent_al_Edit, tline, self.accepted_start_t, self.accepted_end_t,self.accepted_mT, self.m_Flag, currentsugg_mT)
			#self.ssent_O_Edit.setText(out_source)
			#self.tsent_O_Edit.setText(out_target)
		else:
			self.ssent_O_Edit.setText("")# sline) #+" ||| |||")
			self.tsent_O_Edit.setText("") #tline) # +" ||| |||")
			self.ssent_al_Edit.setText("")
			self.tsent_al_Edit.setText("")

	def get_line_with_accepted_term_pairs(self, sent_O_Edit, sent_al_Edit, line, accepted_start, accepted_end, accepted_m, m_Flag, currentsugg):
		global blackColor, greenColor, blueColor, redColor, whitecolor

		sent_O_Edit.setText("")
		sent_al_Edit.setText("")
		for index, word in enumerate(line.split()):
			if index in accepted_start:
				indices = []
				Tlen = {}
				for i, x in enumerate(accepted_start):
					if x == index:
						term_id = i
						term_len = accepted_end[i] - accepted_start[i]
						Tlen[term_len] = term_id

				for key in sorted(Tlen.keys(), reverse=True):
					indices.append(Tlen[key])
				for termID in indices: #[::-1]:
					sent_O_Edit.insertPlainText("<trm"+str(termID+1)+"> ")
				if currentsugg in indices and m_Flag > 0:
					sent_O_Edit.setTextColor(redColor)
					sent_O_Edit.setTextBackgroundColor(greenColor)

			sent_O_Edit.insertPlainText(word+ " ")
			if index in accepted_end:
				indices = [i for i, x in enumerate(accepted_end) if x == index]

				if currentsugg in indices and m_Flag > 0:
					sent_O_Edit.setTextColor(blackColor)
					sent_O_Edit.setTextBackgroundColor(whitecolor)

				for termID in indices:
					sent_O_Edit.insertPlainText("</trm> ")

		sent_O_Edit.insertPlainText("\n")
		#print ("In output:")
		#print (self.accepted_mS)
		#print (accepted_m)
		IDs = ""
		Monosynset= ""
		for i in range (len(accepted_start)):
			term = ' '.join(line.split()[accepted_start[i]:accepted_end[i]+1]).strip()
			#print (term)
			#print (accepted_start[i])
			#print (accepted_end[i])
			IDs += str(i+1)+":"+str(accepted_start[i]) +"-"+ str(accepted_end[i]) +" "
			if (accepted_start[i], accepted_end[i]) in accepted_m:
				accepted_syn = accepted_m[(accepted_start[i], accepted_end[i])]
				synstr = ""
				for syn in accepted_syn:
					synstr+=syn+":"
				Monosynset+= "trm"+str(i+1)+":" + term +":" +synstr.strip(':') + "\n"
			else:
				Monosynset+= "trm"+str(i+1)+":" + term + ":\n"
			

		sent_O_Edit.insertPlainText(Monosynset)
		sent_al_Edit.insertPlainText(IDs)
		
	def get_additional_sugg_from_TermFinder(self):
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		if self.TermFinderTableStoT:
			(tmp_start_s, tmp_end_s, tmp_start_t, tmp_end_t, termsugg_mono) = self.search_s_t_line (sline, tline, self.TermFinderTableStoT, self.max_sterm_len_ST, len(self.TermFinderTableStoT))
			self.get_remaining_sugg((tmp_start_s, tmp_end_s, tmp_start_t, tmp_end_t))
			self.append_synonyms(termsugg_mono, 1)
			#print ("ret 3", termsugg_mono)
			#print("monosugg 3: ",self.termsugg_mono_S, self.termsugg_mono_T)
		if self.TermFinderTableTtoS:
			(tmp_start_t, tmp_end_t, tmp_start_s, tmp_end_s, termsugg_mono) = self.search_s_t_line (tline, sline, self.TermFinderTableTtoS, self.max_sterm_len_TS, len(self.TermFinderTableTtoS))
			self.get_remaining_sugg((tmp_start_s, tmp_end_s, tmp_start_t, tmp_end_t))
			self.append_synonyms(termsugg_mono, 0)
			#print ("ret 4", termsugg_mono)
			#print("monosugg 4: ",self.termsugg_mono_S, self.termsugg_mono_T)

	def append_synonyms(self, termsugg_mono, side):
		for key in termsugg_mono:
			synsugg = termsugg_mono[key]
			if side == 0:
				if key in self.termsugg_mono_S:
					self.termsugg_mono_S[key].extend(synsugg)
				else:
					self.termsugg_mono_S[key] = synsugg
			else:
				if key in self.termsugg_mono_T:
					self.termsugg_mono_T[key].extend(synsugg)
				else:
					self.termsugg_mono_T[key] = synsugg	

	def get_remaining_sugg(self, tup):
		for i in range(len(tup[0])):
			if tup[0][i] not in self.termsugg_start_s or tup[1][i] not in self.termsugg_end_s or tup[2][i] not in self.termsugg_start_t or tup[3][i] not in self.termsugg_end_t:
				self.termsugg_start_s.append(tup[0][i])
				self.termsugg_end_s.append(tup[1][i])
				self.termsugg_start_t.append(tup[2][i])
				self.termsugg_end_t.append(tup[3][i])

	def mark_for_colouring(self, tup):
		#self.from_annot
		for i in range(len(tup[0])):
			fourind = (tup[0][i], tup[1][i], tup[2][i], tup[3][i]) 
			if fourind not in self.from_annot:
				self.from_annot.append(fourind)

	def skip_this_sugg(self):
		self.no_skip_tp+=1
		self.currentsugg+=1
		self.start_annotation()

	def reject_this_sugg_for_future(self):
		self.no_rej_tp+=1
		if self.currentsugg < len(self.termsugg_start_s):
			self.rejected_start_s.append(self.termsugg_start_s[self.currentsugg])
			self.rejected_end_s.append(self.termsugg_end_s[self.currentsugg])
			self.rejected_start_t.append(self.termsugg_start_t[self.currentsugg])
			self.rejected_end_t.append(self.termsugg_end_t[self.currentsugg])

		self.currentsugg+=1
		self.start_annotation()

	def accept_this_sugg(self):
		self.no_acc_tp+=1
		if self.currentsugg < len(self.termsugg_start_s):
			self.accept_term_pair(self.termsugg_start_s[self.currentsugg], self.termsugg_end_s[self.currentsugg], self.termsugg_start_t[self.currentsugg], self.termsugg_end_t[self.currentsugg])
		self.currentsugg+=1
		self.start_annotation()

	def accept_term_pair(self, s_s, s_e, t_s, t_e):
		self.accepted_start_s.append(s_s)
		self.accepted_end_s.append(s_e)
		self.accepted_start_t.append(t_s)
		self.accepted_end_t.append(t_e)

	def skip_syn_sugg(self):
		#print ("Skipping:>>>")
		self.no_skip_m+=1
		self.next_synonym_or_next_term(0)

	def restart_showing_current_term(self):
		self.currentsyn_m_pos = -1

		#print (self.curr_rejected_mS)
		#print (self.accepted_mS)
		#print (self.accepted_mT)
		if self.m_Flag == 1:
			key = (self.accepted_start_s [self.currentsugg_mS], self.accepted_end_s[self.currentsugg_mS])
			if key in self.accepted_mS:
				self.no_acc_m1 -= len(self.accepted_mS [key])
				#del self.accepted_mS [key]
			if key in self.curr_rejected_mS:
				self.no_rej_m1 -= len(self.curr_rejected_mS [key])
				del self.curr_rejected_mS [key]
		if self.m_Flag == 2:
			key = (self.accepted_start_t[self.currentsugg_mT], self.accepted_end_t[self.currentsugg_mT])
			if key in self.accepted_mT:
				self.no_acc_m2 -= len(self.accepted_mT [key])
				#del self.accepted_mT [key]
			if key in self.curr_rejected_mT:
				self.no_rej_m2 -= len(self.curr_rejected_mT [key])
				del self.curr_rejected_mT [key]	

		self.displayOutput()
		self.next_synonym_or_next_term(0)
		#print (self.curr_rejected_mS)
		#print (self.accepted_mS)
		#print (self.accepted_mT)

	def reject_syn_sugg(self):
		synonym = self.textbox_syn.text()
		if synonym != '-':
			term = self.textbox_term.text()
			if self.m_Flag == 1:
				self.no_rej_m1+=1
				#self.curr_rejected_mS [(term, synonym)] = 1
				key = (self.accepted_start_s [self.currentsugg_mS], self.accepted_end_s[self.currentsugg_mS])
				if key in self.curr_rejected_mS:
					self.curr_rejected_mS [key].append(synonym)
				else:
					self.curr_rejected_mS [key] = [synonym]
			if self.m_Flag == 2:
				self.no_rej_m2+=1
				#self.curr_rejected_mT [(term, synonym)] = 1
				key = (self.accepted_start_t[self.currentsugg_mT], self.accepted_end_t[self.currentsugg_mT])
				if key in self.curr_rejected_mT:
					self.curr_rejected_mT[key].append(synonym)
				else:
					self.curr_rejected_mT[key] = [synonym]
			
		self.next_synonym_or_next_term(0)
		#self.displayOutput()

	def next_synonym_or_next_term(self, status):
		#print ("Go to next_synonym_or_next_term>>>")
		#print (self.currentsyn_m_pos)
		#print (self.currentsynlist)
		#print(self.currentsugg_mS)
		#print(self.currentsugg_mT)
		#print ("in: ", status)
		#print (self.currentsynlist)
		#print (self.currentsyn_m_pos)
		syn = self.currentsynlist[self.currentsyn_m_pos]
		print ("here: start: next_synonym_or_next_term")
		if status == 0:
			self.currentsyn_m_pos+=1
		if status == 1 and syn != "-":
			self.currentsyn_m_pos+=1
		
		if self.currentsyn_m_pos < len(self.currentsynlist):
			self.show_syn()
		else:
			if self.m_Flag == 1:
				self.currentsugg_mS+=1
				#if self.currentsugg_mS < len(self.accepted_start_s):
				print ("here: start 2 src: next_synonym_or_next_term")
				self.point_current_sugg_S_or_T()
			elif self.m_Flag == 2:
				self.currentsugg_mT+=1
				#if self.currentsugg_mT < len(self.accepted_start_t):
				print ("here: start 3 src: next_synonym_or_next_term")
				self.point_current_sugg_S_or_T()

	def accept_syn_sugg(self):
		#print ("Storing:>>>")
		synonym = self.textbox_syn.text()
		# use self.currentsugg_mS or self.currentsugg_mT for storing
	
		#skip if syn is '-'	
		if synonym == '-':
			self.next_synonym_or_next_term(0)
		else:
			if self.m_Flag == 1:
				self.no_acc_m1+=1
				if self.currentsugg_mS < len(self.accepted_start_s):
					key = (self.accepted_start_s[self.currentsugg_mS], self.accepted_end_s[self.currentsugg_mS])
					if key in self.accepted_mS:
						self.accepted_mS [key].append(synonym)
					else:
						self.accepted_mS [key] = [synonym]
			elif self.m_Flag == 2:
				self.no_acc_m2+=1
				if self.currentsugg_mT < len(self.accepted_start_t):
					key = (self.accepted_start_t[self.currentsugg_mT], self.accepted_end_t[self.currentsugg_mT])
					if key in self.accepted_mT:
						self.accepted_mT [key].append(synonym)
					else: 
						self.accepted_mT [key] = [synonym]

			#print (self.accepted_mS)
			#print (self.accepted_mT)
			print ("here 1: end os accept_syn_sugg")
			self.displayOutput()
			self.next_synonym_or_next_term(1)

	def show_syn(self):
		term = "END"
		syn = "END"
		if self.m_Flag == 1:
			if self.currentsugg_mS < len(self.accepted_start_s):
				sline = self.bicor[self.sent_ID-1][0]
				term = ' '.join(sline.split()[self.accepted_start_s[self.currentsugg_mS]:self.accepted_end_s[self.currentsugg_mS]+1]).strip()

		if self.m_Flag == 2:
			if self.currentsugg_mT < len(self.accepted_start_t):
				tline = self.bicor[self.sent_ID-1][1]
				term = ' '.join(tline.split()[self.accepted_start_t[self.currentsugg_mT]:self.accepted_end_t[self.currentsugg_mT]+1]).strip()
		if self.currentsyn_m_pos < len(self.currentsynlist):
			syn = self.currentsynlist[self.currentsyn_m_pos]
			#self.textbox_sugg_no.setText(":"+str(self.currentsyn_m_pos))

		#print ("Showing in textbox >>")
		#print ("term: "+term)
		#print ("syn: "+syn)
		setcolor = 0
		if term == "END" and syn == "END":
			setcolor = 0
		else:
			if self.m_Flag == 1:
				if term.lower() in self.mono_syn_S and syn.lower() in self.mono_syn_S[term.lower()]:
					setcolor = 1
				else:
					setcolor = 2
			if self.m_Flag == 2:
				if term.lower() in self.mono_syn_T and syn.lower() in self.mono_syn_T[term.lower()]:
					setcolor = 1
				else:
					setcolor = 2
		if setcolor == 1:
			self.textbox_term.setStyleSheet("color: green;")
			self.textbox_syn.setStyleSheet("color: green;")
			self.textbox_sugg_no.setStyleSheet("color: green;")
		if setcolor == 2:
			self.textbox_term.setStyleSheet("color: blue;")
			self.textbox_syn.setStyleSheet("color: blue;")
			self.textbox_sugg_no.setStyleSheet("color: blue;")
		if setcolor == 0:
			self.textbox_term.setStyleSheet("color: balck;")
			self.textbox_syn.setStyleSheet("color: black;")
			self.textbox_sugg_no.setStyleSheet("color: black;")

		if self.m_Flag == 1:
			self.textbox_sugg_no.setText("s"+str(self.currentsugg_mS+1))
		if self.m_Flag == 2:
			self.textbox_sugg_no.setText("t"+str(self.currentsugg_mT+1))
		if setcolor == 0:
			self.textbox_sugg_no.setText("")
			self.m_Flag = 0

		self.textbox_term.setText(term)
		self.textbox_syn.setText(syn)

		self.displayOutput()	

	def set_synlist(self, synlist):
		self.currentsynlist = []
		self.currentsyn_m_pos = 0
		for synonym in synlist:
			if synonym not in self.currentsynlist:
				self.currentsynlist.append(synonym)
		#print (self.currentsynlist)
		#print (self.currentsyn_m_pos)

	def point_current_sugg_S_or_T(self):
		#print ("In: point_current_sugg_S_or_T>>>")
		#print ("ms:"+str(self.currentsugg_mS))
		#print ("mt:"+str(self.currentsugg_mT))
		if self.currentsugg_mS < len(self.accepted_start_s):
			#print ("sourcesugg",self.termsugg_mono_S)
			#print ("src:",self.accepted_start_s[self.currentsugg_mS], self.accepted_end_s[self.currentsugg_mS])
			if (self.accepted_start_s[self.currentsugg_mS], self.accepted_end_s[self.currentsugg_mS]) in self.termsugg_mono_S:
				synlist = self.termsugg_mono_S[(self.accepted_start_s[self.currentsugg_mS], self.accepted_end_s[self.currentsugg_mS])]
				#if len(synlist) == 0:
				#	self.currentsugg_mS+=1
				#	self.point_current_sugg_S_or_T()
				#else:
				print ("source", synlist)
				self.m_Flag = 1
				self.set_synlist(synlist)
				self.show_syn()
			else:
				self.currentsugg_mS+=1
				self.point_current_sugg_S_or_T()
		elif self.currentsugg_mT < len(self.accepted_start_t):
			#print ("targetsugg",self.termsugg_mono_T)
			#print ("src_sugg",self.accepted_start_t[self.currentsugg_mT], self.accepted_end_t[self.currentsugg_mT])
			if (self.accepted_start_t[self.currentsugg_mT], self.accepted_end_t[self.currentsugg_mT]) in self.termsugg_mono_T:
				synlist = self.termsugg_mono_T[(self.accepted_start_t[self.currentsugg_mT], self.accepted_end_t[self.currentsugg_mT])]
				#if len(synlist) == 0:
				#	self.currentsugg_mT+=1
				#	self.point_current_sugg_S_or_T()
				#else:
				print ("target", synlist)
				self.m_Flag = 2
				self.set_synlist(synlist)
				self.show_syn()
			else:
				self.currentsugg_mT+=1
				self.point_current_sugg_S_or_T()
		else:
			#print ("Here")
			self.currentsynlist = []
			self.show_syn()
			self.start_editing_syn_sugg()

	def amend_addl_syns_from_synfiles(self):
		#self.mono_syn_S
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		for i in range (len(self.accepted_start_s)):
			sterm = ' '.join(sline.split()[self.accepted_start_s[i]:self.accepted_end_s[i]+1]).strip().lower()
			tterm = ' '.join(tline.split()[self.accepted_start_t[i]:self.accepted_end_t[i]+1]).strip().lower()

			#print (sterm)
			#print (tterm+"\n\n")

			#self.accepted_mS [key].append(synonym)

			# Source side 
			key = (self.accepted_start_s[i], self.accepted_end_s[i])
			if sterm in self.mono_syn_S:
				for syn in self.mono_syn_S[sterm]:
					if key in self.accepted_mS:
						self.accepted_mS [key].append(syn)
					else:
						self.accepted_mS [key] = [syn]
			
				if key in self.termsugg_mono_S:
					synlist = []
					for syn in self.termsugg_mono_S [key]:
						if syn not in self.accepted_mS [key]:
							synlist.append(syn)
					self.termsugg_mono_S [key] = synlist
					#self.termsugg_mono_S [ key ] = self.mono_syn_S [ sterm ] + self.termsugg_mono_S [ key ]
				#else:
				#	self.termsugg_mono_S [ key ] = self.mono_syn_S [ sterm ]

			#Target
			key = (self.accepted_start_t[i], self.accepted_end_t[i])
			if tterm in self.mono_syn_T:
				for syn in self.mono_syn_T [tterm]:
					if key in self.accepted_mT:
						self.accepted_mT [key].append(syn)
					else:
						self.accepted_mT [key] = [syn]

				if key in self.termsugg_mono_T:
					synlist = []
					for syn in self.termsugg_mono_T [key]:
						if syn not in self.accepted_mT [ key ]:
							synlist.append(syn)
					self.termsugg_mono_T [key] = synlist
					#self.termsugg_mono_T [ key ] = self.mono_syn_T [ tterm ] + self.termsugg_mono_T[ key ]
				#else:
				#	self.termsugg_mono_T [ key ] = self.mono_syn_T [ tterm ]  
	
		#print ("end == amend_addl_syns_from_synfile")	
	def create_depth_2_dic_S_T(self):
		self.mono_syn_S_depth_2 = {}
		self.mono_syn_T_depth_2 = {}
		
		for key in self.mono_syn_S:
			synset = set(self.mono_syn_S[key])
			sec_set = set()
			for syn in synset:
				sec_set.add(syn)
				sec_set |=  set(self.mono_syn_S[syn])

			sec_set.remove(key)
			self.mono_syn_S_depth_2[key] = sec_set

		for key in self.mono_syn_T:
			synset = set(self.mono_syn_T[key])
			sec_set = set()
			for syn in synset:
				sec_set.add(syn)
				sec_set |=  set(self.mono_syn_T[syn])
			
			sec_set.remove(key)

			self.mono_syn_T_depth_2[key] = sec_set

		print (len(self.mono_syn_S_depth_2))
		print (len(self.mono_syn_T_depth_2))

		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		print ("Before depth2",self.termsugg_mono_S)
		print ("Before depth2",self.termsugg_mono_T)	
		for i in range (len(self.accepted_start_s)):
			sterm = ' '.join(sline.split()[self.accepted_start_s[i]:self.accepted_end_s[i]+1]).strip().lower()
			tterm = ' '.join(tline.split()[self.accepted_start_t[i]:self.accepted_end_t[i]+1]).strip().lower()

			if sterm in self.mono_syn_S_depth_2:
				depth_2_synlist_s = list(self.mono_syn_S_depth_2[sterm])
				key = (self.accepted_start_s[i], self.accepted_end_s[i])
				if key in self.termsugg_mono_S:
					self.termsugg_mono_S[key].extend(list(depth_2_synlist_s))
				else:
					self.termsugg_mono_S[key] = list(depth_2_synlist_s)
			if tterm in self.mono_syn_T_depth_2:
				depth_2_synlist_t = list(self.mono_syn_T_depth_2[tterm])
				key = (self.accepted_start_t[i], self.accepted_end_t[i])
				if key in self.termsugg_mono_T:
					self.termsugg_mono_T[key].extend(list(depth_2_synlist_t))
				else:
					self.termsugg_mono_T[key] = list(depth_2_synlist_t)
		print ("After depth 2:",self.termsugg_mono_S)
		print ("After depth2",self.termsugg_mono_T)

	def remove_dup_mono_sugg(self):
		
		sline = self.bicor[self.sent_ID-1][0]
		tline = self.bicor[self.sent_ID-1][1]

		#print ("before remoning dup:", self.termsugg_mono_S)
		#print ("before remoning dup:", self.termsugg_mono_T)

		new_sugg_dic_src = {}
		new_sugg_dic_trg = {}
		for i in range (len(self.accepted_start_s)):
			sterm = ' '.join(sline.split()[self.accepted_start_s[i]:self.accepted_end_s[i]+1]).strip().lower()
			tterm = ' '.join(tline.split()[self.accepted_start_t[i]:self.accepted_end_t[i]+1]).strip().lower()

			if (self.accepted_start_s[i], self.accepted_end_s[i]) in self.termsugg_mono_S:
				synlist_s = self.termsugg_mono_S[(self.accepted_start_s[i], self.accepted_end_s[i])]
				#print (sterm, synlist_s)
				syn_fl = []
				for synonym in synlist_s:
					if (sterm, synonym.lower()) in self.rejected_mono_s:
						continue
					if synonym.lower() == sterm:
						continue
					if synonym.lower() in syn_fl:
						continue
					syn_fl.append(synonym.lower())
				#syn_fl.append("-")
				if syn_fl != []:
					#del self.termsugg_mono_S[(self.accepted_start_s[i], self.accepted_end_s[i])]
				#else:
					new_sugg_dic_src [(self.accepted_start_s[i], self.accepted_end_s[i])] = syn_fl
				#	self.termsugg_mono_S[(self.accepted_start_s[i], self.accepted_end_s[i])] = syn_fl
			#else:
			#	self.termsugg_mono_S[(self.accepted_start_s[i], self.accepted_end_s[i])] = [] #["-"]

			if (self.accepted_start_t[i], self.accepted_end_t[i]) in self.termsugg_mono_T:
				synlist_t = self.termsugg_mono_T[(self.accepted_start_t[i], self.accepted_end_t[i])]

				syn_fl = []
				for synonym in synlist_t:
					if (tterm, synonym.lower()) in self.rejected_mono_t:
						continue
					if synonym.lower() == tterm:
						continue
					if synonym.lower() in syn_fl:
						continue
					syn_fl.append(synonym.lower())
				#syn_fl.append("-")
				if syn_fl != []:
					#del self.termsugg_mono_T[(self.accepted_start_t[i], self.accepted_end_t[i])]
					new_sugg_dic_trg [(self.accepted_start_t[i], self.accepted_end_t[i])] = syn_fl
				#else:
				#	self.termsugg_mono_T[(self.accepted_start_t[i], self.accepted_end_t[i])] = syn_fl
			#else:
			#	self.termsugg_mono_T[(self.accepted_start_t[i], self.accepted_end_t[i])] = [] #["-"]
	
		#for key, val in self.termsugg_mono_S.items():
		#	if val
			
		self.termsugg_mono_S = new_sugg_dic_src
		self.termsugg_mono_T = new_sugg_dic_trg
		#print ("sugg source: ",self.termsugg_mono_S)
		#print ("sugg taregt: ",self.termsugg_mono_T)


	def init_params(self):
		
		self.currentsugg_mS = 0
		self.currentsugg_mT = 0
		self.currentsyn_m_pos = 0
		self.currentsynlist = []
		self.accepted_mT = {}
		self.accepted_mS = {}
		self.m_Flag = 0
		self.curr_rejected_mS = {}
		self.curr_rejected_mT = {}

		#######################
		self.no_acc_m1 = 0
		self.no_rej_m1 = 0

		self.no_acc_m2 = 0
		self.no_rej_m2 = 0
		self.no_skip_m = 0
	

	def start_syn_sugg(self):
		self.monostarttime = time.time()
		#print ("Starting:")
		self.init_params()

		if self.checkBox.isChecked():
			self.create_depth_2_dic_S_T()

		self.amend_addl_syns_from_synfiles()
		self.remove_dup_mono_sugg()

		#if not self.nosyn_check:
		self.point_current_sugg_S_or_T()
		#else:
		#self.displayOutput_2()

	def select_and_accept_this_sugg(self):
		self.no_new_tp+=1
		self.select_and_accept_sugg()
		self.currentsugg+=1
		self.start_annotation()
	
	def select_and_accept_sugg(self):
		#edited_source = self.ssent_Edit.textCursor().selectedText()  #.toPlainText()
		#edited_target = self.tsent_Edit.textCursor().selectedText()  #.toPlainText()

		(s_s, s_e) = self.identify_start_end_of_selection(self.bicor[self.sent_ID-1][0], self.ssent_Edit.textCursor().selectionStart(), self.ssent_Edit.textCursor().selectionEnd()-1)
		(t_s, t_e) = self.identify_start_end_of_selection(self.bicor[self.sent_ID-1][1], self.tsent_Edit.textCursor().selectionStart(), self.tsent_Edit.textCursor().selectionEnd()-1)

		if s_s == -1 or s_e == -1 or t_s == -1 or t_e == -1:
			print ("Terms Selected Incorrectly >> Error!")
		else:
			#sline = self.bicor[self.sent_ID-1][0]
			#sterm = ' '.join(sline.split()[s_s:s_e+1]).strip().lower()
			#tline = self.bicor[self.sent_ID-1][1]
			#tterm = ' '.join(tline.split()[t_s:t_e+1]).strip().lower()
			
			
			#self.append_synonyms(termsugg_mono, 1)
			self.accept_term_pair(s_s, s_e, t_s, t_e)

		#print ("orig:")#(str(len(self.bicor[self.sent_ID-1][1])))	
		#print (str(self.ssent_Edit.textCursor().selectionStart()))
		#print (str(self.ssent_Edit.textCursor().selectionEnd()))

		#print (str(s_s)+" "+str(s_e)+" "+str(t_s)+" "+str(t_e))

		#print ("["+edited_source+"]")
		#print (edited_target)

	def identify_start_end_of_selection(self, line, si, ei):
		linelist = line.split()
		ret_start = -1
		ret_end = -1
		#print (str(si)+" "+str(ei))
		while line[si] == ' ':
			si+=1
		while line[ei] == ' ':
			ei-=1
		#print (str(si)+" "+str(ei))
		spos=0
		epos=0
		for index, word in enumerate(linelist):
			if spos == si:
				ret_start = index
			spos+=len(word) + 1
			epos+=len(word) - 1
			if epos == ei:
				ret_end = index
			epos+=2

		return (ret_start, ret_end)
			

	def trace_in_target(self, lowerwords_target, tline_tags, tterm):
		#print ("in trace:===>>")
		ids =[]
		tterml = tterm.split()
		for index in range(len(lowerwords_target)):
			if index + len(tterml) <= len (lowerwords_target):
				#print (lowerwords_target[index:index+len(tterml)])
				if self.nottagged(tline_tags[index:index+len(tterml)]) and self.equal(lowerwords_target[index:index+len(tterml)], tterml):
					i = index
					while i < index + len(tterml):
						ids.append(i)
						i+=1
					#print (ids)
					break

		return ids

	def equal(self, A, B):
		for index, word in enumerate(A):
			if word != B[index]:
				return False
		return True

	def nottagged(self, tterm_tags):
		for tag in tterm_tags:
			if tag == '':
				return True
		return False

	def get_suggested_annoatations(self):
                sline = self.bicor[self.sent_ID-1][0]
                tline = self.bicor[self.sent_ID-1][1]
                d = 0
                if self.termTableST:
                    d= len(self.termTableST)
                (self.termsugg_start_s, self.termsugg_end_s, self.termsugg_start_t, self.termsugg_end_t, termsugg_mono) = self.search_s_t_line (sline, tline, self.termTableST, self.max_ann_sterm_len_ST, d) #len(self.termTableST))
                self.append_synonyms(termsugg_mono, 1)
		#print ("ret 1", termsugg_mono)
		#print("monosugg: 1 ",self.termsugg_mono_S, self.termsugg_mono_T)
                self.mark_for_colouring((self.termsugg_start_s, self.termsugg_end_s, self.termsugg_start_t, self.termsugg_end_t))

                d = 0
                if self.termTableTS:
                    d= len(self.termTableTS)
                (tmp_start_t, tmp_end_t, tmp_start_s, tmp_end_s, termsugg_mono) = self.search_s_t_line (tline, sline, self.termTableTS, self.max_ann_sterm_len_TS, d) # len(self.termTableTS))
                self.get_remaining_sugg((tmp_start_s, tmp_end_s, tmp_start_t, tmp_end_t))	
                self.append_synonyms(termsugg_mono, 0)
		#print ("ret 2", termsugg_mono)
		#print("monosugg: 2",self.termsugg_mono_S, self.termsugg_mono_T)
                self.mark_for_colouring((tmp_start_s, tmp_end_s, tmp_start_t, tmp_end_t))

	def search_s_t_line(self, sline, tline, termTable, width_n, termTableSize):
		#print (sline+" <> "+tline)
		#print ("start")

		words_target = tline.split()
		lowerwords_target = tline.lower().split()
		words = sline.split()
		lowerwords = sline.lower().split()

		wordslen = len(lowerwords)
		wordslen_target = len(lowerwords_target)

		#sline_tags = [""] * wordslen
		tline_tags = [""] * wordslen_target

		start = 0

		temp_start_s = []
		temp_end_s = []
		temp_start_t = []
		temp_end_t = []
		termsugg_mono = {}

		while start < wordslen:
			old_start = start
			index = 0
			for i in range(width_n,-1,-1):
				end = start + i
				if end > wordslen:
					continue
				source = u" ".join(lowerwords[start:end]).strip()
				if len(source)+1 > termTableSize:
					continue
				if not source in termTable[len(source)]:
					continue
				#if not self.nottagged(sline_tags[start:end]):
				#	continue
				#print (source)
				notfound = True
				ttermlist = termTable[len(source)][source] #.sort(key=len)
				ttermlist = sorted(ttermlist, key=len, reverse=True)
				#print (ttermlist)
				for tterm in ttermlist: 
					#print(">>>>>>>"+tterm)
					ids = self.trace_in_target(lowerwords_target, tline_tags, tterm)
					if len(ids) > 0:
						temp_start_s.append (start)
						temp_end_s.append (end-1)
						temp_start_t.append (ids[0])
						temp_end_t.append (ids[-1])
						key = (ids[0],ids[-1])
						#print (start,end-1)
						#print (source)
						#print (key)
						#print (tterm)
						#print (ttermlist)
						while tterm in ttermlist:
							ttermlist.remove(tterm)
						#print (ttermlist)
						if len(ttermlist) > 0:
							ttermlist2 = []
							for tsyn in ttermlist:
								if tsyn.lower() not in ttermlist2:
									ttermlist2.append(tsyn.lower())
							if tterm.lower() in ttermlist2:
								ttermlist2.remove(tterm.lower())
							if len(ttermlist2) > 0:
								if key in termsugg_mono:
									termsugg_mono[key].extend(ttermlist2)
								else:
									termsugg_mono[key] = ttermlist2

						notfound = False
						for i in ids:
							tline_tags[i] = "TERM"
						#print ("tsugg")
						#print (termsugg_mono)
						#print (source +" <> " +tterm)
						#print (start, end-1, ids[0], ids[-1])
						break
				if notfound:
					continue
				
				start = end
				break
			if start == old_start:
				start = old_start + 1
		#print ("retstart")
		#print (termsugg_mono)
		#print ("retend")
		return (temp_start_s, temp_end_s, temp_start_t, temp_end_t, termsugg_mono)


	def initUI(self):
		global bicor

		ssent  = QLabel('Source Text', self)
		tsent  = QLabel('Target Text', self)
		ssent.move(250, 10)
		tsent.move(790, 10)

		self.ssent_Edit = QTextEdit(self) 
		self.tsent_Edit = QTextEdit(self)
		self.ssent_Edit.resize(500,130)
		self.tsent_Edit.resize(500,130)
		self.ssent_Edit.move(50, 40)
		self.tsent_Edit.move(600, 40)

		self.ssent_Edit.setReadOnly(True)
		self.tsent_Edit.setReadOnly(True)

		self.loadb = QPushButton('Load sentence ID:', self)
		self.loadb.move(50, 180)
		self.loadb.resize(180, 25)
		self.loadb.setStyleSheet("background-color: red")
		
		self.textbox = QLineEdit(self)
		self.textbox.move(250, 180)
		self.textbox.resize(60, 25)

		self.nextb = QPushButton('Skip this pair for now!', self)
		self.nextb.move(600, 180)
		self.nextb.resize(180, 25)

		self.reject_sent_b = QPushButton('Reject this pair!', self)
		self.reject_sent_b.move(900, 180)
		self.reject_sent_b.resize(180, 25)
		self.reject_sent_b.setStyleSheet("background-color: red")

		##############################################################

		sugg_label1  = QLabel('Term Suggestions:', self)
		sugg_label1.move(220, 245)

		self.textbox_sugg1 = QLineEdit(self)
		self.textbox_sugg1.move(370, 240)
		self.textbox_sugg1.resize(60, 25)

		self.skipb = QPushButton('Skip', self)
		self.skipb.move(460, 240)

		self.rejectb = QPushButton('Reject for Future', self)
		self.rejectb.move(460, 275)
		self.rejectb.resize(130, 25)

		self.acceptb = QPushButton('Accept', self)
		self.acceptb.move(550, 240)

		self.selectb = QPushButton("Accept Selected Terms", self)
		self.selectb.move(650, 240)
		self.selectb.resize(180, 25)

		####################################################################
		ssentO  = QLabel('Annotated Source Text', self)
		tsentO  = QLabel('Annotated Target Text', self)
		ssentO.move(200, 300)
		tsentO.move(700, 300)

		#h=220
		#if not self.nosyn_check:
		#h=150	
		self.ssent_O_Edit = QTextEdit(self)
		self.tsent_O_Edit = QTextEdit(self)
		self.ssent_O_Edit.resize(500, 175)
		self.tsent_O_Edit.resize(500, 175)
		self.ssent_O_Edit.move(50, 330)
		self.tsent_O_Edit.move(600, 330)

		self.ssent_al_Edit = QTextEdit(self)
		self.tsent_al_Edit = QTextEdit(self)
		self.ssent_al_Edit.resize(500, 30)
		self.tsent_al_Edit.resize(500, 30)
		self.ssent_al_Edit.move(50, 520)
		self.tsent_al_Edit.move(600, 520)


		self.save_annot = QPushButton("SAVE", self)
		self.save_annot.move(1080, 630)
		#self.save_annot.resize(70, 40)
		self.save_annot.setStyleSheet("background-color: red")

		######################################################
		#if not self.nosyn_check:
		mono  = QLabel('Term:', self)
		mono.move(50, 565) #495)

		self.textbox_sugg_no = QLineEdit(self)
		self.textbox_sugg_no.move(100, 565) #490)
		self.textbox_sugg_no.resize(35, 25)
		
		self.restart_showing = QPushButton('Reset', self)
		self.restart_showing.move(143, 565) #490)
		self.restart_showing.resize(55, 25)
	
		mono  = QLabel('Suggestion from TermFinder:', self)
		mono.move(280, 565) #495)
		
		self.textbox_term = QLineEdit(self)
		self.textbox_term.move(50, 595) #520)
		self.textbox_term.resize(200, 25)

		self.textbox_syn = QLineEdit(self)
		self.textbox_syn.move(280, 595) #520)
		self.textbox_syn.resize(200, 25)

		self.skipsyn = QPushButton('Skip', self)
		self.skipsyn.move(280, 630) #590) #560)
		
		self.rejectsyn = QPushButton('Reject for Future', self)
		self.rejectsyn.move(140, 630 ) #600)
		self.rejectsyn.resize(130,25)

		self.acceptsyn = QPushButton('Accept', self)
		self.acceptsyn.move(380, 630) # 590) #560)


		self.startsyn = QPushButton('Start', self)
		self.startsyn.move(45, 630) # 590)  #560)
			#self.startsyn.resize(40,50)
		self.startsyn.setStyleSheet("background-color: red")

		self.checkBox = QCheckBox('Syn+', self)
		self.checkBox.move(45, 660)
		self.checkBox.setChecked(True)

		#self.editsyn = QPushButton('StartSynEdit', self)
		#self.editsyn.move(540, 630) # 505)
		#self.startsyn.resize(40,50)
		#self.editsyn.setStyleSheet("background-color: red")

		##########################################################
		self.elim_dup = QPushButton('Sort Terminology and \n Synonym Dictionaries', self)
		self.elim_dup.move(900, 230)
		self.elim_dup.resize(170,50)
		self.elim_dup.setStyleSheet("background-color: blue")

		self.elim_dup.clicked.connect(self.eliminate_duplicates_from_term_bases)

		########################
		self.loadb.clicked.connect(self.getSentenceID)
		self.nextb.clicked.connect(self.displayNextSentence)
		self.reject_sent_b.clicked.connect(self.RejectThisSent_displayNextSentence)

		self.skipb.clicked.connect(self.skip_this_sugg)
		self.acceptb.clicked.connect(self.accept_this_sugg)
		self.selectb.clicked.connect(self.select_and_accept_this_sugg)
		self.rejectb.clicked.connect(self.reject_this_sugg_for_future)


		self.startsyn.clicked.connect(self.start_syn_sugg)
		#self.editsyn.clicked.connect(self.start_editing_syn_sugg)
		#if not self.nosyn_check:
		self.save_annot.clicked.connect(self.save_this_annotated_sent_pairs)
		self.restart_showing.clicked.connect(self.restart_showing_current_term)
		self.skipsyn.clicked.connect(self.skip_syn_sugg)
		self.acceptsyn.clicked.connect(self.accept_syn_sugg)
		self.rejectsyn.clicked.connect(self.reject_syn_sugg)

		# Load available sentence starting from ID := 1
		self.get_next_available_ID()
		self.loadSentence()

		self.setGeometry(20, 20, 1200, 690)
		self.setWindowTitle('TermMarker')
		self.show()
        
if __name__ == '__main__':

	argparser = argparse.ArgumentParser()
	argparser.add_argument('-sf', '--sfile', required=True, help='Source file')
	argparser.add_argument('-tf', '--tfile', required=True, help='Target file')
	argparser.add_argument('-termstf', '--termstfile', required=False, help='Terminology (bilingual) s t1 t2 t3 t4')
	argparser.add_argument('-termtsf', '--termtsfile', required=False, help='Terminology (bilingual) t s1 s2 s3 s4')
	#argparser.add_argument('--nosyncheck', dest='nosyn_check',  action="store_true", required=False, help='Second Mode')
		
	args = argparser.parse_args()

	app = QApplication(sys.argv)
	ta = TermAnnotator(args.sfile, args.tfile, args.termstfile, args.termtsfile) #, args.nosyn_check)
	sys.exit(app.exec_())
