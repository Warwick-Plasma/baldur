
Clone the repo using:
git clone git@cfsa-pmw.warwick.ac.uk:SDF/baldur.git

add to your .bashrc :
export PYTHONPATH="${PYTHONPATH}:/home/duncan/baldur"
export PYTHONSTARTUP=~/.pythonrc 

create a file called .pythonrc with this inside:
# modules
import sdf_helper as sh
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import numpy as np
import glob
import sys, os
import tkinter as tk
from tkinter import ttk
import odin_plot as op
import import_sdf as isdf
import menu as mu
print('Modules and constants loaded')

you will need SDF installed to proceed any further:
cd SDF/C
make
cd ../FORTRAN
make
cd ../utilities
./build

Now in start ipython with:
ipython

and open the options menu with:
Python 3.6.3 |Anaconda, Inc.| (default, Oct 13 2017, 12:02:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 6.1.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: import menu as mu

In [2]: mu.options(use_analysis=True)

In order for the code to run: 
-plotting must be active, you must have an sdf file (*.sdf) in your current working directory, and you need an up-to-date sdf_helper. 
When running with Odin certain dump masks are required to set use_analysis = True.

To save videos ffmpeg must be installed, try:
sudo apt  install ffmpeg
