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
VERSION = "v1.0.1"
FILE = "2zc1.fasta"
PH = 7
HYDROLIMIT = 0
NUMCUT = 1
GRAPH1 = 1
GRAPH2 = 2
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
STYLE = "seaborn-muted"

#amino acid class
class aminoAcid:
	def __init__(self, name, hphob2, hphob7):
		self.name =  name
		self.hphob2 = hphob2
		self.hphob7= hphob7
#instantiating all amino acid objects
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
	try:
		file = open(filename, 'r')
	except:
		errorMessage("Filename not found!", 180)
		return
	sequenceString = file.read().replace(" ","").replace("\n", "")
	file.close()
	sequenceList = list(sequenceString)
	remove = False
	removeRest1 = False
	removeRest2 = False
	codecount = 0
	for char in sequenceList[:]:
		print(char)
		if removeRest2:
			sequenceList.remove(char)
		elif char == '>':
			if removeRest1:
				removeRest2 = True
			removeRest1 = True
			sequenceList.remove(char)
			remove = True
		elif removeRest1 == False:
			break
		elif remove:
			if codecount == 0:
				if char == 'N':
					codecount = 1
			elif codecount == 1:
				if char == 'C':
					codecount = 2
			elif codecount == 2:
				if char == 'E':
					remove = False
					print("false!")
			else:
				codecount = 0
			sequenceList.remove(char)
		#print("remove: ", remove)
		#print("removeRest1: ", removeRest1)
		#print("removeRest2: ", removeRest2)
		#print("codecount: ", codecount)

	file = open(filename, 'w')
	file.write(''.join(sequenceList))
	file.close()
	#print(sequenceList)
	return sequenceList
class Application(ttk.Frame):
	def __init__(self, master = None):
		ttk.Frame.__init__(self, master)
		self.pack()
		self.initialize()
	def initialize(self):
		self.graphExists = False
		self.genFrame = Tk.Frame(master = root)
		self.genFrame.pack(side = Tk.LEFT)
		self.inputsFrame = ttk.Frame(master = self.genFrame)
		self.inputsFrame.pack(side = Tk.TOP, padx= 0, pady = 6)
		self.fileFrame = Tk.Frame(master = self.inputsFrame, bg = "#bbc3cc")
		self.fileFrame.pack(side = Tk.TOP, pady = 0, padx = 0)
		self.fileLabel = Tk.Label(master = self.fileFrame, text = "File Name:              ", bg = "#bbc3cc")
		self.fileLabel.pack(side = Tk.LEFT, padx = 3, pady = 3)
		self.fileTkVar = Tk.StringVar()
		self.fileEntry = ttk.Entry(master = self.fileFrame, width = 10, textvariable = self.fileTkVar)
		self.fileEntry.pack(side = Tk.LEFT, padx = 4, pady = 5)
		self.fileTkVar.set(FILE)
		self.numCutFrame = Tk.Frame(master = self.inputsFrame, bg = "#9fbfdf")
		self.numCutFrame.pack(side = Tk.TOP, pady = 0)
		self.numCutLabel = Tk.Label(master = self.numCutFrame, text = "Number of Cutoffs:        ", bg = "#9fbfdf")
		self.numCutLabel.pack(side = Tk.LEFT, padx = 4, pady = 3)
		self.numCutTkVar = Tk.IntVar()
		self.numCutCombobox = ttk.Combobox(master = self.numCutFrame, values = [1,2,3,4], width = 2, textvariable = self.numCutTkVar,
			 state = "readonly")
		self.numCutCombobox.bind("<<ComboboxSelected>>", lambda e: self.updateCutoff(self))
		self.numCutTkVar.set(NUMCUT)
		self.currNumCutoffs = 1
		self.numCutCombobox.pack(side = Tk.LEFT, padx = 4, pady = 3)
		self.pHFrame = Tk.Frame(master = self.inputsFrame, bg = "#bbc3cc")
		self.pHFrame.pack(side = Tk.TOP, pady = 0) 
		self.pHLabel = Tk.Label(master = self.pHFrame, text = "pH:                                    ", bg = "#bbc3cc")
		self.pHLabel.pack(side = Tk.LEFT, padx = 3, pady = 3)
		self.pHTkVar = Tk.IntVar()
		self.pHCombobox = ttk.Combobox(master = self.pHFrame, values = [2,7], textvariable = self.pHTkVar, state = "readonly", width = 2)
		self.pHFrame = Tk.Frame(master = self.inputsFrame, bg = "#bbc3cc")
		self.pHTkVar.set(PH)
		self.pHCombobox.pack(side = Tk.LEFT, padx = 5, pady = 3)
		self.graph1Frame = Tk.Frame(master = self.inputsFrame, bg = "#9fbfdf")
		self.graph1Frame.pack(side = Tk.TOP, pady = 0) 
		self.graph1Label = Tk.Label(master = self.graph1Frame, text = "Graph 1:                           ", bg = "#9fbfdf#")
		self.graph1Label.pack(side = Tk.LEFT, padx = 4, pady = 3)
		self.graph1TkVar = Tk.IntVar()
		self.graph1Combobox = ttk.Combobox(master = self.graph1Frame, values = [1,2,3,4,5], textvariable = self.graph1TkVar, state = "readonly", width = 2)
		self.graph1TkVar.set(GRAPH1)
		self.graph1Combobox.pack(side = Tk.LEFT, padx = 5, pady = 3)
		self.graph2Frame = Tk.Frame(master = self.inputsFrame, bg = "#bbc3cc")
		self.graph2Frame.pack(side = Tk.TOP, pady = 0) 
		self.graph2Label = Tk.Label(master = self.graph2Frame, text = "Graph 2:                           ", bg = "#bbc3cc")
		self.graph2Label.pack(side = Tk.LEFT, padx = 4, pady = 3)
		self.graph2TkVar = Tk.IntVar()
		self.graph2Combobox = ttk.Combobox(master = self.graph2Frame, values = [1,2,3,4,5], textvariable = self.graph2TkVar, state = "readonly", width = 2)
		self.graph2TkVar.set(GRAPH2)
		self.graph2Combobox.pack(side = Tk.LEFT, padx = 5, pady = 3)
		self.cutFrameList = []
		self.cutFrameList.append(self.numCutFrame)
		self.hydroLimitFrame = Tk.Frame(master = self.inputsFrame, bg = "#9fbfdf")
		self.hydroLimitFrame.pack(side = Tk.TOP, pady = 0)
		self.hydroLimitLabel = Tk.Label(master = self.hydroLimitFrame, text = "Hydrophobicity Cutoff 1:", bg = "#9fbfdf")
		self.hydroLimitLabel.pack(side = Tk.LEFT, padx = 0, pady = 3)
		self.hydroLimitTkVarList = []
		hydroLimitTkVar = Tk.IntVar()
		self.hydroLimitTkVarList.append(hydroLimitTkVar)
		self.hydroLimitEntry = ttk.Entry(master = self.hydroLimitFrame, width = 5, textvariable = hydroLimitTkVar)
		hydroLimitTkVar.set(HYDROLIMIT)
		self.hydroLimitEntry.pack(side = Tk.LEFT, padx = 5, pady = 5)
		self.currColor = "blue"
		self.buttonFrame = Tk.Frame(master = self.genFrame)
		self.buttonFrame.pack(side = Tk.TOP)
		self.analyzeButton = ttk.Button(master = self.buttonFrame, text = "Analyze", width = 8, command = lambda: self.analyze())
		self.analyzeButton.pack(side = Tk.TOP, pady = 6)
		self.canvasFrame = Tk.Frame(master = root)
		self.canvasFrame.pack(side = Tk.LEFT, expand = 1)
		self.startCanvas = Tk.Canvas(master = self.canvasFrame, width = 550, height = 330, bd = 0, relief = "ridge", highlightthickness = 0)
		self.startCanvas.pack(side = Tk.TOP)
		self.startCanvas.configure(bg = "#dde8f1")
		self.startCanvas.create_text(280, 330 / 3 + 30, text = "Welcome to ProteinAnalyzer! You must have a .txt file with your desired sequence\n " + 
			"in the same directory as the executable. To analyze, type in your .txt filename, \n" + " adjust the inputs, and press Analyze.")
		#self.photo = Tk.PhotoImage(file = "gfp.png")
		#self.startCanvas.create_image(546 / 2, 326 / 2, image = self.photo)
	def updateCutoff(event, self):
		newNumCuts = int(self.numCutTkVar.get())
		print("currNumCuts:", self.currNumCutoffs)
		print("newNumCuts: ", newNumCuts)
		
		if newNumCuts < self.currNumCutoffs:
			frameID = self.currNumCutoffs
			while frameID != newNumCuts:
				self.cutFrameList[frameID - 1].destroy()
				del self.cutFrameList[-1]
				if self.currColor == "blue":
					self.currColor = "gray"
				else:
					self.currColor = "blue"
					print("Removed ", frameID)
					print("\n")
				frameID -= 1
		elif newNumCuts > self.currNumCutoffs:
			for frameID in range(self.currNumCutoffs + 1, newNumCuts + 1):
				if self.currColor == "blue":
					hydroLimitFrame = Tk.Frame(master = self.inputsFrame, bg = "#bbc3cc")
					hydroLimitFrame.pack(side = Tk.TOP, pady = 0)
					hydroLimitLabel = Tk.Label(master = hydroLimitFrame, text = "Hydrophobicity Cutoff %i:" %frameID, bg = "#bbc3cc")
					hydroLimitLabel.pack(side = Tk.LEFT, padx = 0, pady = 3)
					hydroLimitTkVar = Tk.IntVar()
					self.hydroLimitTkVarList.append(hydroLimitTkVar)
					hydroLimitEntry = ttk.Entry(master = hydroLimitFrame, width = 5, textvariable = hydroLimitTkVar)
					hydroLimitTkVar.set(HYDROLIMIT)
					hydroLimitEntry.pack(side = Tk.LEFT, padx = 5, pady = 5)
					self.cutFrameList.append(hydroLimitFrame)
					self.currColor = "gray"
				else:
					hydroLimitFrame = Tk.Frame(master = self.inputsFrame, bg = "#9fbfdf")
					hydroLimitFrame.pack(side = Tk.TOP, pady = 0)
					hydroLimitLabel = Tk.Label(master = hydroLimitFrame, text = "Hydrophobicity Cutoff %i:" %frameID, bg = "#9fbfdf")
					hydroLimitLabel.pack(side = Tk.LEFT, padx = 0, pady = 3)
					hydroLimitTkVar = Tk.IntVar()
					self.hydroLimitTkVarList.append(hydroLimitTkVar)
					hydroLimitEntry = ttk.Entry(master = hydroLimitFrame, width = 5, textvariable = hydroLimitTkVar)
					hydroLimitTkVar.set(HYDROLIMIT)
					hydroLimitEntry.pack(side = Tk.LEFT, padx = 5, pady = 5)
					self.cutFrameList.append(hydroLimitFrame)
					self.currColor = "blue"
				print("\n")
		else:
			return
			print("\n")
			#assert False, "newNumcuts equals currNumCuts, should not happen!"
		self.currNumCutoffs = newNumCuts

	def analyze(self):
		filename = self.fileTkVar.get() + ".txt"
		#print(filename)
		self.sequenceList = parseFile(filename)
		#print(sequenceList)
		self.sortedAminoList = self.sortAminoAcids()
		self.plotDistributions()
	def plotDistributions(self):
		style.use(STYLE)
		font = {'fontname':'Helvetica'}
		if self.graphExists:
			self.subplot1.clear()
			self.subplot2.clear()
			#self.canvas.get_tk_widget().destroy()
			#self.toolbar.destroy()
		#print(self.sortedAminoList)
		acidID1 = int(self.graph1TkVar.get()) - 1
		acidID2 = int(self.graph2TkVar.get()) - 1
		histData1 = self.getHistogramData(acidID1)
		histData2 = self.getHistogramData(acidID2)
		#print(histData2)
		if not self.graphExists	:
			self.startCanvas.destroy()
			self.plotFigure = Figure(figsize=(5.5, 3.3), dpi=100, facecolor = "#dde8f1")
			self.subplot1 = self.plotFigure.add_subplot(121)
			self.subplot2 = self.plotFigure.add_subplot(122)
			self.subplot1.set_color_cycle(COLORARRAY)
			self.subplot2.set_color_cycle(COLORARRAY)
			self.subplot1.tick_params(labelsize = 7)
			self.subplot2.tick_params(labelsize = 7)
		binwidth = 1
		maximum = max(max(histData1), max(histData2))
		x1, y1, _ = self.subplot1.hist(histData1, color = COLORARRAY[0], normed = True, 
			bins=range(min(histData1), maximum + binwidth, binwidth))
		self.subplot1.set_ylabel("Normalized Frequency", labelpad=5, fontsize = 8, **font)
		self.subplot1.set_xlabel("Hydrophilic Block Size" , labelpad = 5, fontsize = 8, **font)
		ylim1 = self.subplot1.get_ylim()
		self.subplot2.hist(histData2, color = COLORARRAY[1], normed = True, 
			bins=range(min(histData2), maximum + binwidth, binwidth))
		self.subplot2.set_ylabel("Normalized Frequency", labelpad=5, fontsize = 8, **font)
		self.subplot2.set_xlabel("Hydrophobic Block Size" , labelpad = 5, fontsize = 8, **font)
		ylim2 = self.subplot2.get_ylim()
		maxLim = max(max(ylim1), max(ylim2))
		print("ylim1: ", ylim1)
		print("ylim2: ", ylim2)
		self.subplot1.set_ylim([0, maxLim])
		self.subplot2.set_ylim([0, maxLim])
		self.plotFigure.tight_layout()
		# A tk.DrawingArea
		#imbedding matplotlib graph onto canvas
		if not self.graphExists:
			self.canvas = FigureCanvasTkAgg(self.plotFigure, master = self.canvasFrame)
		else:
			self.canvas.draw()
		self.canvas.show()
		
		#Imbedding matplotlib toolbar onto canvas
		if not self.graphExists:
			#self.toolbar = NavigationToolbar2TkAgg(self.canvas, root)
			#self.toolbar.update()
			self.canvas._tkcanvas.pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 1)
			self.canvas.get_tk_widget().configure(bg = "black", bd = 0, selectbackground = "black")
			self.canvas.get_tk_widget().pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 1)
		self.graphExists = True
		print("height: ", self.canvasFrame.winfo_height())
		print("width: ", self.canvasFrame.winfo_width())
	def unpackTkVarList(self, tkVarList):
		intList = []
		for tkVar in tkVarList:
			intList.append(int(tkVar.get()))
		return intList
	def sortAminoAcids(self):
		sortedList = []
		self.numCuts = int(self.numCutTkVar.get())
		self.pH = int(self.pHTkVar.get())
		self.limitList = self.unpackTkVarList(self.hydroLimitTkVarList)
		self.limitList = self.limitList[0:self.numCuts]
		self.limitList = sorted(set(self.limitList))
		self.numCuts = len(self.limitList)
		print("limitList: ", self.limitList)
		for aminoAcid in self.sequenceList:
			aminoAcidObj = getAminoAcid(aminoAcid)
			added = False
			limitID = 0
			for limit in self.limitList:
				if not added:
					if self.pH == 2:
						if aminoAcidObj.hphob2 <= limit:
							sortedList.append(limitID)
							added = True
					if self.pH == 7:
						if aminoAcidObj.hphob7 <= limit:
							sortedList.append(limitID)
							added = True
				limitID += 1
			if not added:
				sortedList.append(limitID)
		print(sortedList)

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
def errorMessage(message, width):
	#Toplevel parameters
	top = Tk.Toplevel()
	top.grab_set()
	top.wm_title("Error")
	top.geometry("%dx%d%+d%+d" % (width, 70, 250, 125))
	#Message
	msg = Tk.Message(master = top, text = message, width = 500)
	msg.pack(side = Tk.TOP, pady = 5)
	#OK button to exit
	exitButton = ttk.Button(master = top, text = "Ok", command = top.destroy, width = 7)
	exitButton.pack(side = Tk.TOP, pady = 5)
root = Tk.Tk()
root.style = ttk.Style()
root.style.theme_use("vista")
root.wm_title("Protein Analyzer %s" % VERSION)
app = Application(master = root)
app.mainloop()
