""" Inputs a single-component waveform (trace or stream not certain yet) and 
	outputs the times of the signals.

	Pick objects might be overkill for this project, but BAM will use them, 
	so I've put them in

	Author: Luke McFadden
"""
import numpy as np

class Pick:

	def __init__(self, first_pt, last_pt):
		""" Make a pick object defined as the outer limits of the region

		Example Usage:
		A = Pick(1.0, 2.0)
		print(A)
		>>> Pick object defined from 1.00 - 2.00 s

		Arguments:
		first_pt [float] - the initial time of the pick
		last_pt [float] - the ending time of the pick

		"""

		try:
			first_pt = float(first_pt)
			last_pt = float(last_pt)
		except ValueError:
			print("[ERROR] Pick must be defined with two floats!")
			first_pt = np.nan
			last_pt = np.nan

		if first_pt > last_pt:
			first_pt, last_pt = last_pt, first_pt

		self.first_pt = first_pt
		self.last_pt = last_pt


	def __str__(self):

		""" Prints the pick boundaries for testing and display
		"""

		return "Pick object defined from {:.2f} - {:.2f} s".format(self.first_pt, self.last_pt)

	def __add__(self, other):
		'''
		Returns a pick object which envelopes the outer maxima of both picks

		Example Usage:
		A = Pick(1, 2)
		B = Pick(2, 3)

		A + B -> Pick(1, 3) 	

		Arguments:
		self [Pick Obj] - the first pick object to be added to the range
		other [Pick Obj] - the second pick object to be added to the range

		Returns [Pick Obj] - the first_pt is the minimum of the two first_pt inputs and the
							 last_pt is the maximum of the two last_pt inputs
		'''

		first_pt = min(self.first_pt, other.first_pt)
		last_pt = max(self.last_pt, other.last_pt)

		return Pick(first_pt, last_pt)

	def __sub__(self, other):

		''' Returns the difference between the two first_pt of the inputs

		Example Usage:
		A = Pick(1, 2)
		B = Pick(2, 3)

		B - A -> 1

		Arguments:
		self [Pick Obj] - the first pick object to be subtracted
		other [Pick Obj] - the second pick object to be subtracted

		Returns [float] - the difference between the first_pts
		'''

		return self.first_pt - other.first_pt

	def getLen(self):
		'''
		Returns the length of the pick (difference between first and last point)

		Example Usage:

		A = Pick(1, 3)
		A.getLen()
		>>> 2

		'''

		return self.last_pt - self.first_pt
		

	

if __name__ == "__main__":

	pass

	A = Pick(1, 2)
	B = Pick(2, 3)

	print("For testing: Error message -> pass, crash -> fail")
	C = Pick("A", "B")

	D = A + B
	E = A - B

	assert (D.first_pt == min(A.first_pt, B.first_pt))
	assert (D.last_pt == max(A.last_pt, B.last_pt))
	assert (E == A.first_pt - B.first_pt)
	assert np.isnan(C.first_pt)
	assert np.isnan(C.last_pt)

	F = Pick(2, 1)
	assert F.first_pt == 1 and F.last_pt == 2

	assert A.getLen() == 1
	assert D.getLen() == 2

	print("Pick testing done!")