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


def diff(i, j): # I think this is where I will decide if the two match or not
	if i == j:
		val = -3  # they're the same, they match
	else:
		val = 1  # they're different so we have to sub
	return val


class GeneSequencing:

	def __init__( self ):
		pass
	
# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	#TODO if we have time @ end - time complexities
	def align( self, seq1, seq2, banded, align_length):   #seq1 and 2 are just the strings we want to run the prog on
		self.banded = banded
		self.MaxCharactersToAlign = align_length   #how many of the characters we want to use
		sub1 = seq1[:align_length]
		sub2 = seq2[:align_length]

		#TODO how do I make a list of lists that holds ints
		self.matrixVals = [[]] # list of lists containing the values
		self.matrixPrev = [[]] # list of lists containing the prev node
		# E = {x,y : (val, prev)}  TA suggested not doing a dictionary but to do 2 lists

 		#self.edit(sub1, sub2)   # e is the matrix i want
# dummy code
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 		score = random.random()*100;
# 		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq1), align_length, ',BANDED' if banded else '')
# 		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################
		score = 1  #the last element of the matrixVals

		#TODO figure out how to recourse out of the matrix and make the first 100 of each string to return
		alignment1 = 1
		alignment2 = 1
		
		return {'align_cost' : score, 'seqi_first100' : alignment1, 'seqj_first100' : alignment2}

	def edit(self, string1, string2):
		if self.banded:
			bandWidth = MAXINDELS
		else:
			bandWidth = self.MaxCharactersToAlign

		# E will be our matrix I make
		for i in range(bandWidth):  # a bunch of inserts along first row
			E[i, 0] = i * INDEL   #TODO how do i index and update 2d matricies
			#TODO I also need to update the prev pointers for these- it would be just the prev x val
		for j in range(bandWidth):  # more inserts along column
			E[0, j] = j * INDEL
		for i in range(1, string1):
			for j in range(1, string2):
				E[i,j] = min(diff(i,j) + E[i-1, j-1], 5 + E[i, j-1], 5 + E[i-1, j]) # (diagonal(match, or sub): -3 or 1 , down(indel): 5, right(indel): 5)
			# TODO accessing elements, i need to update the vals - are these the right vals/changes?
			# TODO how do I have "min" return the one it picked so that I can update the prev vals matrix as well


		# return E[[i], [j]] I dont think I need to return anything since im storing the matricies in the class
