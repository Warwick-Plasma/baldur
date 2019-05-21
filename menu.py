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
from functools import partial



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
  
  fig = plt.figure(num=1,figsize=(10,8),facecolor='white')
  ax1 = plt.axes()
  setattr(ax1, 'cbar', 'None')
  move_figure(fig, 550, 10)

  def callbackFunc(event):
    grid_colour = grid_variable.get()
    use_polar = polar_variable.get()
    sdf_num = slider1.get()
    var_name = combo1.get()
    
    dat = isdf.use_sdf(sdf_num, pathname, use_analysis = use_analysis, istart = istart)
    op.snapshot(dat, ax1, var_name = var_name, grid_colour = grid_colour,
                      use_polar = use_polar)
  
  filename1 = 'test.pdf'
  def save_pdf():
    fig.savefig(filename1)

  app = tk.Tk() 
  app.geometry('450x200+10+10')

  # slider - time
  label_slider1 = tk.Label(app, text = "Select sdf number:")
  label_slider1.grid(column=0, row=0)

  slider1 = tk.Scale(app, from_ = istart, to = RunCounter, tickinterval=100, orient=tk.HORIZONTAL, command=callbackFunc, length  = 300, resolution = 1.0)
  slider1.grid(column=1, row=0)
  slider1.set(23)

  # Combo box - variable
  labelTop_combo1 = tk.Label(app, text = "Select variable:")
  labelTop_combo1.grid(column=0, row=1)

  combo1 = ttk.Combobox(app, values = dat.variables)
  combo1.grid(column=1, row=1)
  combo1.current(3)
  
  # check box - grid
  grid_variable = tk.StringVar(app)
  grid_button = tk.Checkbutton(app, text="grid", variable = grid_variable, onvalue="black", offvalue="None")
  grid_button.deselect()
  grid_button.grid(column=0, row=2)
  
  # check box - polar coordinates
  polar_variable = tk.BooleanVar(app)
  polar_button = tk.Checkbutton(app, text="polar coordinates", variable = polar_variable, onvalue=True, offvalue=False)
  polar_button.deselect()
  polar_button.grid(column=0, row=3)
  
  # button - save fig as pdf
  print_button = tk.Button(app, text="Save .pdf", command = save_pdf)
  print_button.grid(column=1, row=2)
  
  # slider key bindings
  def leftKey(event):
    sdf_num = slider1.get()
    slider1.set(sdf_num - 1)

  def rightKey(event):
    sdf_num = slider1.get()
    slider1.set(sdf_num + 1)
  
  app.bind('<Left>', leftKey)
  app.bind('<Right>', rightKey)
  combo1.bind("<<ComboboxSelected>>", callbackFunc)
  grid_button.bind("<ButtonRelease-1>", callbackFunc)
  polar_button.bind("<ButtonRelease-1>", callbackFunc)

  app.mainloop()



def main():
	"""
	"""


if __name__ == "__main__":
	main()
