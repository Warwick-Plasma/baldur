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
  """Move figure's upper left corner to pixel (x, y)
  https://yagisanatode.com/2018/02/23/how-do-i-change-the-size-and-position-of-the-main-window-in-tkinter-and-python-3/
  """
  backend = matplotlib.get_backend()
  if backend == 'TkAgg':
    f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
  elif backend == 'WXAgg':
    f.canvas.manager.window.SetPosition((x, y))
  else:
    # This works for QT and GTK
    # You can also use window.setGeometry
    f.canvas.manager.window.move(x, y)



def options(*args, **kwargs):
  use_analysis = kwargs.get('use_analysis', False)
  root = tk.Tk()
  my_gui = menu_GUI(root, use_analysis)
  root.mainloop()



class menu_GUI:
  def __init__(self, app, use_analysis):
    self.app = app
    self.use_analysis = use_analysis
    app.title("A simple GUI") 
    app.geometry('450x200+10+10')
    
    plt.ion()
    plt.close('all')
    
    self.pathname = os.path.abspath(os.getcwd())
    runs = glob.glob1(self.pathname,"*.sdf")
    RunCounter = len(runs)
    run_array = np.zeros(RunCounter)
    for ir in range(0, RunCounter-1):
      run_name = runs[ir]
      run_num = int(run_name[:-4])
      run_array[ir] = run_num
    run_array = sorted(run_array)
    self.istart = int(run_array[0])
    self.iend = int(run_array[-1])
    sdf_num = int(run_array[0])

    dat = isdf.use_sdf(sdf_num, self.pathname, use_analysis = self.use_analysis, istart = self.istart)
    self.fig = plt.figure(num=1, figsize=(6,6), facecolor='white')
    move_figure(self.fig, 700, 10)
    self.ax1 = plt.axes()
    setattr(self.ax1, 'cbar', 'None')

    cs = 1
    self.dat0 = isdf.get_data_all(dat, self.istart, self.iend, self.pathname, self.use_analysis, cs)
    self.fig2 = plt.figure(num=2, figsize=(6,6), facecolor='white')
    move_figure(self.fig2, 700, 1400)
    self.ax2 = plt.axes()
    self.ax3 = op.empty_lineout(self.dat0, self.fig2, self.ax2)
    
    self.reset_grid_variable = tk.BooleanVar(app)
    self.reset_grid_variable.set(True)
    
    # slider - time
    self.label_slider1 = tk.Label(app, text = "Select sdf number:")
    self.label_slider1.grid(column=0, row=0)

    self.slider1 = tk.Scale(app, from_ = self.istart, to = self.iend, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.callbackFunc,
                            length  = 300, resolution = 1.0)
    self.slider1.grid(column=1, row=0)
    self.slider1.set(self.istart)

    # Combo box - variable
    self.labelTop_combo1 = tk.Label(app, text = "Select variable:")
    self.labelTop_combo1.grid(column=0, row=1)

    self.combo1 = ttk.Combobox(app, values = dat.variables)
    self.combo1.grid(column=1, row=1)
    self.combo1.current(3)
  
    # check box - grid
    self.grid_variable = tk.StringVar(app)
    self.grid_button = tk.Checkbutton(app, text="grid", variable = self.grid_variable,
                                      onvalue="black", offvalue="None")
    self.grid_button.deselect()
    self.grid_button.grid(column=0, row=2)
  
    # check box - polar coordinates
    self.polar_variable = tk.BooleanVar(app)
    self.polar_button = tk.Checkbutton(app, text="polar coordinates",
                                       variable = self.polar_variable,
                                       onvalue=True, offvalue=False)
    self.polar_button.deselect()
    self.polar_button.grid(column=0, row=3)
    
    # check box - anisotropies
    self.anisotropies_variable = tk.BooleanVar(app)
    self.anisotropies_button = tk.Checkbutton(app, text="View anisotropies",
                                       variable = self.anisotropies_variable,
                                       onvalue=True, offvalue=False)
    self.anisotropies_button.deselect()
    self.anisotropies_button.grid(column=0, row=4)
  
    # button - save fig as pdf
    self.print_button = tk.Button(app, text="Save .pdf", command = self.save_pdf)
    self.print_button.grid(column=1, row=2)
    
    # button - reset button
    self.reset_button = tk.Button(app, text="Reset zoom")
    self.reset_button.grid(column=1, row=3)
  
    self.app.bind('<Left>', self.leftKey)
    self.app.bind('<Right>', self.rightKey)
    self.combo1.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.grid_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.polar_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.anisotropies_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.reset_button.bind("<Button-1>", self.callbackFunc1)
    
  def callbackFunc1(self, event):
    self.reset_grid_variable.set(True)
    self.callbackFunc(event)

  def callbackFunc(self, event):
    grid_colour = self.grid_variable.get()
    use_polar = self.polar_variable.get()
    sdf_num = self.slider1.get()
    var_name = self.combo1.get()
    reset_axis = self.reset_grid_variable.get()
    view_anisotropies = self.anisotropies_variable.get()

    dat = isdf.use_sdf(sdf_num, self.pathname, use_analysis = self.use_analysis, istart = self.istart)
    op.snapshot(dat, self.fig, self.ax1, var_name = var_name,
        grid_colour = grid_colour, use_polar = use_polar,
        reset_axis = reset_axis, view_anisotropies = view_anisotropies)
    
    op.lineout(self.dat0, self.fig2, self.ax2, self.ax3, var_name, sdf_num-self.istart)

    self.reset_grid_variable.set(False)

  def save_pdf(self):
    filename1 = 'test.pdf'
    self.fig.savefig(filename1)
  
  def leftKey(self, event):
    sdf_num = self.slider1.get()
    self.slider1.set(sdf_num - 1)

  def rightKey(self, event):
    sdf_num = self.slider1.get()
    self.slider1.set(sdf_num + 1)



def main():
	"""
	"""


if __name__ == "__main__":
	main()
