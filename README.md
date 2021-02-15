#Pre-requisite
Clone and build SDF repo
git clone --recursive git@cfsa-pmw.warwick.ac.uk:SDF/SDF.git
navigate to utilities and run .build

Anaconda is useful but python and many modules are required for SDF and Baldur

#Install
Clone the repo using:
git clone git@cfsa-pmw.warwick.ac.uk:SDF/baldur.git

Add to your .bashrc :
export PYTHONPATH="${PYTHONPATH}:/home/duncan/baldur"

#Run
To run use:
python [install location]/baldur/menu.py 0 0
from the location of an Odin or Epoch output file ie and .sdf file. The zeroes at the end can be changed to ones to turn on analysis and time history plot

#Optional

For use with ipython follow instructions below.

Add to your .bashrc :
export PYTHONSTARTUP=~/.pythonrc

create a file called .pythonrc with this inside:

import sdf_helper as sh
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import odin_plot as op
import import_sdf as isdf
import menu as mu
print('Modules and constants loaded')

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
sudo apt install ffmpeg
