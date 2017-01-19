"""
    Created by Vincent Wu on 1/18/17:
    This program analyzes the hydrophobic and hydrophillic blocks of a protein amino acid sequence
"""
#All neccesary imports
import matplotlib
import random
import math
import json
import os.path
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import style
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from Tkinter import ttk
    import Tkinter.filedialog
    import Tkinter.messagebox
    print(sys.version_info[0], "hi")
else:
    import tkinter as Tk
    from tkinter import ttk
    from tkinter.colorchooser import *
    import tkinter.filedialog
    import tkinter.messagebox
    print(sys.version_info[0])
#global variables
VERSION = "v1.0"
FILE = "5g5y.fasta"
PH = 7
HYDROLIMIT = 0
STYLE = "bmh"
COLOR1 = '#4D4D4D'
COLOR2 = '#5DA5DA'
COLOR3 = '#F15854'
COLOR4 = '#DECF3F'
COLOR5 = '#60BD68'
COLOR6 = '#F17CB0'
COLOR7 = '#B276B2'
COLOR8 = '#FAA43A'
COLORARRAY = [COLOR1, COLOR2, COLOR3, COLOR4, COLOR5, COLOR6, COLOR7, COLOR8]
#amino acid class
class aminoAcid:
	def __init__(self, name, hphob2, hphob7):
		self.name =  name
		self.hphob2 = hphob2
		self.hphob7= hphob7
#instantiating all amino acid onjects
ALA = aminoAcid("A", 47, 41)
ARG = aminoAcid("R", -26, -14)
ASN = aminoAcid("N", -41, -28)
ASP = aminoAcid("D", -18, -55)
CYS = aminoAcid("C",  52, 49)
GLU = aminoAcid("E", 8, -31)
GLN = aminoAcid("Q", -18, -10)
GLY = aminoAcid("G", 0, 0)
HIS = aminoAcid("H", -42, 8)
HYP = aminoAcid("O", None, None)
ILE = aminoAcid("I", 100, 99)
LEU = aminoAcid("L", 100, 97)
LYS = aminoAcid("K", -37, -23)
MET = aminoAcid("M", 74, 74)
PHE = aminoAcid("F", 92, 100)
PRO = aminoAcid("P", -46, -46)
GLP = aminoAcid("U", None, None)
SER = aminoAcid("S", -7, -5)
THR = aminoAcid("T", 13, 13)
TRP = aminoAcid("W", 84, 97)
TYR = aminoAcid("Y", 49, 63)
VAL = aminoAcid("V", 79, 76)
aminoAcidList = [ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, HYP, ILE, LEU, LYS, 
MET, PHE, PRO, GLP, SER, THR, TRP, TYR, VAL]
#given the letter name, returns the correct amino acid object
def getAminoAcid(name):
	for acid in aminoAcidList:
		if acid.name == name:
			return acid
	assert False, "Amino Acid %s not found!" %name
def parseFile(filename):
	file = open(filename, 'r')
	sequenceString = file.read().replace(" ","").replace("\n", "")
	sequenceList = list(sequenceString)
	#print(sequenceList)
	return sequenceList
class Application(ttk.Frame):
	def __init__(self, master = None):
		ttk.Frame.__init__(self, master)
		self.pack()
		self.initialize()
	def initialize(self):
		self.graphExists = False
		self.inputsFrame = ttk.Frame(master = root)
		self.inputsFrame.pack(side = Tk.LEFT, padx= 6, pady = 6)
		self.fileFrame = ttk.Frame(master = self.inputsFrame)
		self.fileFrame.pack(side = Tk.TOP, pady = 3)
		self.fileLabel = ttk.Label(master = self.fileFrame, text = "File Name:              ")
		self.fileLabel.pack(side = Tk.LEFT)
		self.fileTkVar = Tk.StringVar()
		self.fileEntry = ttk.Entry(master = self.fileFrame, width = 10, textvariable = self.fileTkVar)
		self.fileEntry.pack(side = Tk.LEFT)
		self.fileTkVar.set(FILE)
		self.hydroLimitFrame = ttk.Frame(master = self.inputsFrame)
		self.hydroLimitFrame.pack(side = Tk.TOP, pady = 3)
		self.hydroLimitLabel = ttk.Label(master = self.hydroLimitFrame, text = "Hydrophobicity Cutoff: ")
		self.hydroLimitLabel.pack(side = Tk.LEFT)
		self.hydroLimitTkVar = Tk.IntVar()
		self.hydroLimitEntry = ttk.Entry(master = self.hydroLimitFrame, width = 5, textvariable = self.hydroLimitTkVar)
		self.hydroLimitTkVar.set(HYDROLIMIT)
		self.hydroLimitEntry.pack(side = Tk.LEFT)
		self.pHFrame = ttk.Frame(master = self.inputsFrame)
		self.pHFrame.pack(side = Tk.TOP, pady = 3) 
		self.pHLabel = ttk.Label(master = self.pHFrame, text = "pH:                                    ")
		self.pHLabel.pack(side = Tk.LEFT)
		self.pHTkVar = Tk.IntVar()
		self.pHCombobox = ttk.Combobox(master = self.pHFrame, values = [2,7], textvariable = self.pHTkVar, state = "readonly", width = 2)
		self.pHTkVar.set(PH)
		self.pHCombobox.pack(side = Tk.LEFT)
		self.analyzeButton = ttk.Button(master = self.inputsFrame, text = "Analyze", width = 8, command = lambda: self.analyze())
		self.analyzeButton.pack(side = Tk.TOP, pady = 3)
	def analyze(self):
		filename = self.fileTkVar.get() + ".txt"
		#print(filename)
		self.sequenceList = parseFile(filename)
		#print(sequenceList)
		self.sortedAminoList = self.sortAminoAcids()
		self.plotDistributions()
	def plotDistributions(self):
		if self.graphExists:
			self.canvas.get_tk_widget().destroy()
			self.toolbar.destroy()
		self.graphExists = True
		#print(self.sortedAminoList)
		histData1 = self.getHistogramData(0)
		histData2 = self.getHistogramData(1)
		#print(histData2)
		self.plotFigure = Figure(figsize=(5.5, 3.3), dpi=100)
		self.subplot1 = self.plotFigure.add_subplot(121)
		self.subplot2 = self.plotFigure.add_subplot(122)
		self.subplot1.set_color_cycle(COLORARRAY)
		self.subplot2.set_color_cycle(COLORARRAY)
		self.subplot1.tick_params(labelsize = 7)
		self.subplot2.tick_params(labelsize = 7)
		self.subplot1.hist(histData1, color = COLORARRAY[0])
		self.subplot1.set_ylabel("Normalized Frequency", labelpad=5, fontsize = 9)
		self.subplot1.set_xlabel("Hydrophilic Block Size" , labelpad = 0, fontsize = 9)
		self.subplot2.hist(histData2, color = COLORARRAY[1])
		self.subplot2.set_ylabel("Normalized Frequency", labelpad=5, fontsize = 9)
		self.subplot2.set_xlabel("Hydrophobic Block Size" , labelpad = 0, fontsize = 9)
		# A tk.DrawingArea
		#imbedding matplotlib graph onto canvas
		self.canvas = FigureCanvasTkAgg(self.plotFigure, master = root)
		self.canvas.show()
		#Imbedding matplotlib toolbar onto canvas
		self.toolbar = NavigationToolbar2TkAgg(self.canvas, root)
		self.toolbar.update()
		self.canvas._tkcanvas.pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 1)
		self.canvas.get_tk_widget().pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 1)
	def sortAminoAcids(self):
		sortedList = []
		self.pH = int(self.pHTkVar.get())
		self.hydroLimit = int(self.hydroLimitTkVar.get())
		for aminoAcid in self.sequenceList:
			aminoAcidObj = getAminoAcid(aminoAcid)
			if self.pH == 2:
				if aminoAcidObj.hphob2 <= self.hydroLimit:
					sortedList.append(0)
				else:
					sortedList.append(1)
			if self.pH == 7:
				if aminoAcidObj.hphob7 <= self.hydroLimit:
					sortedList.append(0)
				else:
					sortedList.append(1)
		return sortedList
		#returns an array of numbers for each consecutive monomer; to be used in histogram plotting
	def getHistogramData(self, acidID):
		#initializing array to be returned
		histogramData = []
		numConsecutive = 0
		#count through all polymers
		for aminoAcid in self.sortedAminoList: 
			
			#if monomer is not consecutive and monomer before was, add consecutive number to data and reset values
			if aminoAcid != acidID and numConsecutive > 0:
				count = 0
				while count < numConsecutive:
					histogramData.append(numConsecutive)
					count += 1
				numConsecutive = 0
				continue
			#increment consecutive counter by 1 if monomer is consecutive
			if  aminoAcid == acidID:
				numConsecutive += 1
				continue
		return histogramData

root = Tk.Tk()
root.wm_title("Protein Analyzer %s" % VERSION)
app = Application(master = root)
app.mainloop()