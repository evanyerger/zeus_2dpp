import numpy as np
import matplotlib.pylab as plt

def readFile(filename):
	allData  = np.loadtxt(filename)
	plt.imshow(allData,origin='lower')
	plt.colorbar()
	plt.title("Transport")
	plt.show()
readFile("output.txt")


