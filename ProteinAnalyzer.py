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
VERSION = "v1.0"
class Application(ttk.Frame):
	def __init__(self, master = None):
		ttk.Frame.__init__(self, master)
		self.pack()
		self.initialize()
root = Tk.Tk()
root.wm_title("Protein Analyzer %s" % VERSION)
app = Application(master = root)
app.mainloop()
