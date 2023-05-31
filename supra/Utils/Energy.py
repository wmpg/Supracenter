

class EnergyObj:
	def __init__(self):

		self.source_type = "fragmentation"

	def __str__(self):

		A = "Energy Object @ {:} km \n".format(self.height/1000)

		B = "\tSource: {:} \n".format(self.source_type)

		try:
			C = "\tYield: {:.2E} J".format(self.chem_pres)
		except:
			C = "\tYield: {:.2E} J/m".format(self.linear_E)

		return A + B + C
