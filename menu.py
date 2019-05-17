# Duncan Barlow, Odin project, Warwick University, 01/19
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



def move_figure(f, x, y): # cxrodgers
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)



def options():

  plt.ion()
  plt.close('all')
  
  istart = 0
  sdf_num = 10
    
  pathname = os.path.abspath(os.getcwd())
  runs = glob.glob1(pathname,"*.sdf")
  RunCounter = len(runs)-1

  use_analysis = True
  dat = isdf.use_sdf(sdf_num, pathname, use_analysis = use_analysis, istart = istart)
  
  fig = plt.figure(num=1,figsize=(12,8),facecolor='white')
  ax1 = plt.axes()
  setattr(ax1, 'cbar', 'None')
  move_figure(fig, 500, 10)

  def callbackFunc(event):
    sdf_num = slider1.get()
    var_name = combo1.get()
    dat = isdf.use_sdf(sdf_num, pathname, use_analysis = use_analysis, istart = istart)
    op.snapshot(dat, ax1, var_name = var_name)

  app = tk.Tk() 
  app.geometry('400x200+10+10')

  label_slider1 = tk.Label(app, text = "Select sdf number:")
  label_slider1.grid(column=0, row=0)

  slider1 = tk.Scale(app, from_ = istart, to = RunCounter, tickinterval=100, orient=tk.HORIZONTAL, command=callbackFunc)
  slider1.grid(column=1, row=0)
  slider1.set(23)

  labelTop_combo1 = tk.Label(app, text = "Select variable:")
  labelTop_combo1.grid(column=0, row=1)

  combo1 = ttk.Combobox(app,
      values = dat.variables)
  print(dict(combo1))
  combo1.grid(column=1, row=1)
  combo1.current(3)
  
  op.set_axis_lim(dat, ax1, combo1.get())

  combo1.bind("<<ComboboxSelected>>", callbackFunc)

  app.mainloop()



def main():
	"""
	"""


if __name__ == "__main__":
	main()
