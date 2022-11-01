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
	def align( self, seq1, seq2, banded, align_length):   # seq1 and 2 are just the strings we want to run the prog on
		self.banded = banded
		self.MaxCharactersToAlign = align_length   # this is how many of the characters we want to use

		sub1 = seq1[:align_length]
		sub2 = seq2[:align_length]

		#TODO how do I make a list of lists that holds ints
		self.matrixVals = []  # list of lists containing the values
		self.matrixPrev = []  # list of lists containing the prev node

		for i in range(align_length):
			self.matrixVals.append([])
			self.matrixPrev.append([])
			for j in range(align_length):
				self.matrixVals[i].append(None)
				self.matrixPrev[i].append(None)
		self.matrixVals[0][0] = 0

 		self.edit(sub1, sub2)

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


		for i in range(1, bandWidth):  # a bunch of inserts along first row
			self.matrixVals[i][0] = i * INDEL
			self.matrixPrev[i][0] = self.matrixPrev[i-1][0]
		for j in range(1, bandWidth):  # more inserts along column
			self.matrixVals[0][j] = j * INDEL
			self.matrixPrev[0][j] = self.matrixPrev[0][j-1]
		for i in range(1, string1):
			for j in range(1, string2):
				#self.matrixVals[i][j] = min(diff(i,j) + self.matrixVals[i-1][j-1], INDEL + self.matrixVals[i][j-1], INDEL + self.matrixVals[i-1][j]) # (diagonal(match, or sub): -3 or 1 , down(indel): 5, right(indel): 5)
				left = INDEL + self.matrixVals[i-1][j]
				top = INDEL + self.matrixVals[i][j-1]
				diag = diff(i, j) + self.matrixVals[i-1][j-1]

				smallestVal = left
				location = (i-1, j)
				if top < smallestVal:
					smallestVal = top
					location = (i, j-1)
				if diag < smallestVal:
					smallestVal = diag
					location = (i-1, j-1)

				self.matrixVals[i][j] = smallestVal
				self.matrixPrev[i][j] = location
		# return E[[i], [j]] I dont think I need to return anything since im storing the matricies in the class
