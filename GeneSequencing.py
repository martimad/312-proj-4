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


# comparison - O(1)
def diff(i, j): # I think this is where I will decide if the two match or not
	if i == j:
		val = MATCH  # they're the same, they match
	else:
		val = SUB  # they're different so we have to sub
	return val


class GeneSequencing:

	def __init__( self ):
		pass

	#O(nm) TIME AND SPACE unbanded, O(kn) TIME AND SPACE banded
	def align( self, seq1, seq2, banded, align_length):   # seq1 and 2 are just the strings we want to run the prog on

		self.banded = banded
		self.MaxCharactersToAlign = align_length   # this is how many of the characters we want to use

		sub1 = " " + seq1[:align_length]
		sub2 = " " + seq2[:align_length]

		self.matrixVals = []  # list of lists containing the values
		self.matrixPrev = []  # list of lists containing the prev node

		#O(1) - assignments and comparisions
		if banded and abs(len(sub1) - len(sub2)) > 200:
			score = math.inf
			alignment1 = "No Alignment Possible"
			alignment2 = "No Alignment Possible"

		else:
			#O(nm) SPACE - n X m matrix - set up for both
			for i in range(len(sub1)):
				self.matrixVals.append([])
				self.matrixPrev.append([])
				for j in range(len(sub2)):
					self.matrixVals[i].append(None)
					self.matrixPrev[i].append((0, 0))
			self.matrixVals[0][0] = 0

			#O(nm) TIME - explained in func
			self.edit(sub1, sub2)

			#O(1) - assignments
			score = self.matrixVals[-1][-1]  #the last element of the matrixVals

			#O(n) TIME, O(1) space
			alignment1, alignment2 = self.returnStrings(sub1, sub2)

			#O(1) - assignments
			if len(alignment1) > 100:
				alignment1 = alignment1[0:101]
			if len(alignment2) > 100:
				alignment2 = alignment2[0:101]

		return {'align_cost' : score, 'seqi_first100' : alignment1, 'seqj_first100' : alignment2}


	#O(nm) SPACE AND TIME unbanded, O(kn) TIME banded
	def edit(self, string1, string2):

		#O(1) - just assignments, setting how long the interation lengths will be
		if self.banded:
			bandWidth1 = MAXINDELS
			bandWidth2 = MAXINDELS
		else:
			bandWidth1 = len(string1)
			bandWidth2 = len(string2)

		#O(n) goes through each number
		for i in range(bandWidth1):  # a bunch of inserts along first row
			self.matrixVals[i][0] = i * INDEL
			self.matrixPrev[i][0] = (i-1, 0)

		#O(n) goes through each number
		for j in range(bandWidth2):  # more inserts along column
			self.matrixVals[0][j] = j * INDEL
			self.matrixPrev[0][j] = (0, j-1)
		self.matrixPrev[0][0] = None

		#O(kn) time-  banded algorithm
		if self.banded:
			#O(n) - eventually goes through every index in i
			for i in range(1, len(string1)):
				smallestVal = math.inf
				location = (0,0)

				#O(k) - up to k values - depending if its at the beginning or end determines the depth of the reach for the j loop
				if i >= 3:
					for j in range(i-MAXINDELS, min(MAXINDELS + i, len(string2))):
						#left
						if self.matrixVals[i - 1][j] is not None:
							left = INDEL + self.matrixVals[i - 1][j]
							smallestVal = left
							location = (i - 1, j)

						#top
						if self.matrixVals[i][j - 1] is not None:
							top = INDEL + self.matrixVals[i][j - 1]
							if top < smallestVal:
								smallestVal = top
								location = (i, j-1)

						if self.matrixVals[i - 1][j - 1] is not None:
							diag = diff(string1[i], string2[j]) + self.matrixVals[i - 1][j - 1]
							if diag < smallestVal:
								smallestVal = diag
								location = (i - 1, j - 1)

						# smallestVal, location = self.innerLoop(i, j, string1, string2)
						self.matrixVals[i][j] = smallestVal
						self.matrixPrev[i][j] = location

				else:
					for j in range(1, min(MAXINDELS + i, len(string2))):
						#smallestVal, location = self.innerLoop(i, j, string1, string2)
						# left
						if self.matrixVals[i - 1][j] is not None:
							left = INDEL + self.matrixVals[i - 1][j]
							smallestVal = left
							location = (i - 1, j)

						# top
						if self.matrixVals[i][j - 1] is not None:
							top = INDEL + self.matrixVals[i][j - 1]
							if top < smallestVal:
								smallestVal = top
								location = (i, j - 1)

						diag = diff(string1[i], string2[j]) + self.matrixVals[i - 1][j - 1]
						if diag < smallestVal:
							smallestVal = diag
							location = (i - 1, j - 1)

						self.matrixVals[i][j] = smallestVal
						self.matrixPrev[i][j] = location

		#O(nm) unbanded algorithm
		else:
			#O(n) goes through each in i
			for i in range(1, bandWidth1):
				#O(m) goes through each in j
				for j in range(1, bandWidth2):
					left = INDEL + self.matrixVals[i - 1][j]
					top = INDEL + self.matrixVals[i][j - 1]
					diag = diff(string1[i], string2[j]) + self.matrixVals[i - 1][j - 1]

					smallestVal = left
					location = (i - 1, j)

					if top < smallestVal:
						smallestVal = top
						location = (i, j - 1)
					if diag < smallestVal:
						smallestVal = diag
						location = (i - 1, j - 1)

					self.matrixVals[i][j] = smallestVal
					self.matrixPrev[i][j] = location

	#O(n) TIME, O(1) SPACE
	def returnStrings(self, seq1, seq2):
		alignment1 = ""
		alignment2 = ""

		# O(1) - accessing stuff
		coord = (len(seq1) - 1, len(seq2) - 1)  # the very last one in prev array
		seq1ind = len(seq1) - 1  # I just changed this, may need to switch back to self.MaxCharactersToAlign - 1
		seq2ind = len(seq2) - 1

		# O(n) - TIME - iterating through at most i+j items
		while True:
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

		#O(1) - some quick comparisions to get the very last letter
		if not seq1ind == 1:
			alignment1 = seq1[seq1ind] + alignment1
		else:
			alignment1 = '-' + alignment1

		if not seq2ind == 1:
			alignment2 = seq2[seq2ind] + alignment2
		else:
			alignment2 = '-' + alignment2

		return alignment1, alignment2
