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

		sub1 = " " + seq1[:align_length]
		sub2 = " " + seq2[:align_length]

		self.matrixVals = []  # list of lists containing the values
		self.matrixPrev = []  # list of lists containing the prev node

		for i in range(len(sub1)):
			self.matrixVals.append([])
			self.matrixPrev.append([])
			for j in range(len(sub2)):
				self.matrixVals[i].append(None)
				self.matrixPrev[i].append((0, 0))
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
		score = self.matrixVals[-1][-1]  #the last element of the matrixVals

		alignment1, alignment2 = self.returnStrings(sub1, sub2)

		if len(alignment1) > 100:
			alignment1 = alignment1[0:101]
		if len(alignment2) > 100:
			alignment2 = alignment2[0:101]

		return {'align_cost' : score, 'seqi_first100' : alignment1, 'seqj_first100' : alignment2}

	def edit(self, string1, string2):
		if self.banded:
			bandWidth1 = MAXINDELS
			bandWidth2 = MAXINDELS
		else:
			bandWidth1 = len(string1)
			bandWidth2 = len(string2)
		for i in range(bandWidth1):  # a bunch of inserts along first row
			self.matrixVals[i][0] = i * INDEL
			self.matrixPrev[i][0] = (i-1, 0)
		for j in range(bandWidth2):  # more inserts along column
			self.matrixVals[0][j] = j * INDEL
			self.matrixPrev[0][j] = (0, j-1)
		self.matrixPrev[0][0] = None
		for i in range(1, len(string1)):
			for j in range(1, bandWidth2):
				left = INDEL + self.matrixVals[i-1][j]
				top = INDEL + self.matrixVals[i][j-1]
				diag = diff(string1[i], string2[j]) + self.matrixVals[i-1][j-1]

				# then its on the top of the band
				if top is None:
					smallestVal = left
					location = (i-1, j)

				#in center or bottom
				else:
					smallestVal = top
					location = (i, j - 1)

				# if it's in the center
				if left is not None:
					if left < smallestVal:
						smallestVal = left
						location = (i-1, j)

				#all need to compare diag
				if diag < smallestVal:
					smallestVal = diag
					location = (i-1, j-1)

				self.matrixVals[i][j] = smallestVal
				self.matrixPrev[i][j] = location

	def returnStrings(self, seq1, seq2):
		alignment1 = ""
		alignment2 = ""
		coord = (len(seq1) - 1, len(seq2) - 1)  # the very last one in prev array
		seq1ind = len(seq1) - 1  # I just changed this, may need to switch back to self.MaxCharactersToAlign - 1
		seq2ind = len(seq2) - 1

		while True:     # TODO can't get the first character
			prevCoor = self.matrixPrev[coord[0]][coord[1]]
			if prevCoor is None:
				break
			xDist = coord[0] - prevCoor[0]
			yDist = coord[1] - prevCoor[1]

			# if from the diag -- sub or match
			if xDist == 1 and yDist == 1:
				alignment1 = seq1[seq1ind] + alignment1
				seq1ind -= 1
				alignment2 = seq2[seq2ind] + alignment2
				seq2ind -= 1


			# if from the left -- added char to seq1
			if xDist == 1 and yDist == 0:
				alignment2 = '-' + alignment2
				alignment1 = seq1[seq1ind] + alignment1
				seq1ind -= 1

			# if from the top -- added char to seq2
			if xDist == 0 and yDist == 1:
				alignment1 = '-' + alignment1
				alignment2 = seq2[seq2ind] + alignment2
				seq2ind -= 1
			coord = prevCoor

		# TODO hardcoding the last index?
		if not seq1ind == 1:
			alignment1 = seq1[seq1ind] + alignment1
		else:
			alignment1 = '-' + alignment1

		if not seq2ind == 1:
			alignment2 = seq2[seq2ind] + alignment2
		else:
			alignment2 = '-' + alignment2

		return alignment1, alignment2
