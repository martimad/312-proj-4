#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
# elif PYQT_VER == 'PYQT4':
# 	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):   #seq1 and 2 are just the strings we want to run the prog on
		self.banded = banded
		self.MaxCharactersToAlign = align_length   #how many of the characters we want to use
		sub1 = seq1[:align_length]
		sub2 = seq2[:align_length]

		#TODO how do I make a list of lists that holds ints
		self.matrixVals = [[]] # list of lists containing the values
		self.matrixPrev = [[]] # list of lists containing the prev node
		# E = {x,y : (val, prev)}  TA suggested not doing a dictionary but to do 2 lists

 		self.edit(sub1, sub2)   # e is the matrix i want
# dummy code
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 		score = random.random()*100;
# 		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq1), align_length, ',BANDED' if banded else '')
# 		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################
		score = 1
		alignment1 = 1
		alignment2 = 1
		
		return {'align_cost' : score, 'seqi_first100' : alignment1, 'seqj_first100' : alignment2}

	def edit(self, string1, string2):
		# E will be our matrix I make
		for i in range(string1):
			E[i, 0] = i
		for j in range(string2):
			E[0, j] = j
		for i in range(1,string1):
			for j in range(1,string2):
				E[i,j] = min(diff(i,j) + E[i-1, j-1], 1 + E[i, j-1], 1 + E[i-1, j])
		return E[[i], [j]]
