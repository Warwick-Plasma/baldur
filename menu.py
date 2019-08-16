# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.animation import FuncAnimation
import matplotlib.animation as ani
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
  op.check_analysis(use_analysis)
  root = tk.Tk()
  my_gui = snapshot_GUI(root, use_analysis)
  root.mainloop()
  root = tk.Tk()
  my_gui = time_history_GUI(root, use_analysis)
  root.mainloop()



class time_history_GUI:
  def __init__(self, app, use_analysis):
    self.app = app
    self.use_analysis = use_analysis
    app.title("Time history GUI")
    app.geometry('450x200+10+10')
    
    plt.ion()
    plt.close('all')
    
    self.cs = 1
    
    # find sdf files and count
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
    
    # initial data import, needed for variable selection combo box
    dat = isdf.use_sdf(sdf_num, self.pathname, use_analysis = self.use_analysis, istart = self.istart)
    self.dat = isdf.get_data_all(dat, self.istart, self.iend, self.pathname, self.use_analysis, self.cs)
    
    # create empty figures
    self.fig = plt.figure(num=1, figsize=(6,6), facecolor='white')
    move_figure(self.fig, 700, 10)
    self.ax1 = plt.axes()
    setattr(self.ax1, 'cbar', 'None')

    self.fig2 = plt.figure(num=2, figsize=(6,6), facecolor='white')
    move_figure(self.fig2, 10, 1000)
    self.ax2 = plt.axes()
    self.ax3 = op.empty_lineout(self.fig2, self.ax2)
    
    # Combo box - variable
    self.labelTop_combo1 = tk.Label(app, text = "Select variable:")
    self.labelTop_combo1.grid(column=0, row=1)

    self.combo1 = ttk.Combobox(app, values = dat.variables)
    self.combo1.grid(column=1, row=1)
    self.combo1.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.combo1.current(0)
    
    # Combo box - time variable
    self.labelTop_combo2 = tk.Label(app, text = "Select time variable:")
    self.labelTop_combo2.grid(column=0, row=2)

    self.combo2 = ttk.Combobox(app, values = dat.variables_time)
    self.combo2.grid(column=1, row=2)
    self.combo2.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.combo2.current(0)
    
    # Combo box - grid choice
    grid_list = ['default', 'initial', 'cell number']
    self.labelTop_combo3 = tk.Label(app, text = "Select a grid:")
    self.labelTop_combo3.grid(column=0, row=3)

    self.combo3 = ttk.Combobox(app, values = grid_list)
    self.combo3.grid(column=1, row=3)
    self.combo3.bind("<<ComboboxSelected>>", self.callbackFunc1)
    self.combo3.current(0)
    
    # slider - scale up colour
    self.label_slider1 = tk.Label(app, text = "Minimum value on colourbar:")
    self.label_slider1.grid(column=0, row=4)

    self.slider1 = tk.Scale(app, from_ = -10, to = 0, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.callbackFunc,
                            length  = 200, resolution = 0.01)
    self.slider1.grid(column=1, row=4)
    self.slider1.set(0)
    
    # button - reset button
    self.reset_button = tk.Button(app, text="Reset zoom")
    self.reset_button.grid(column=0, row=5)
    self.reset_axis_variable = tk.BooleanVar(app)
    self.reset_axis_variable.set(True)
    
    self.reset_button.bind("<Button-1>", self.callbackFunc1)
  
  def callbackFunc1(self, event):
    self.reset_axis_variable.set(True)
    self.callbackFunc(event)
  
  def callbackFunc(self, event):
    var_name = self.combo1.get()
    var_name2 = self.combo2.get()
    grid_choice = self.combo3.get()
    cbar_upscale = self.slider1.get()
    reset_axis = self.reset_axis_variable.get()
      
    op.time_history(self.dat, self.fig, self.ax1, var_name = var_name, cbar_upscale = cbar_upscale, grid = grid_choice, reset_axis = reset_axis)
      
    op.time_history_lineout(self.dat, self.fig2, self.ax2, self.ax3, var_name = var_name2,  use_analysis = self.use_analysis)

    self.reset_axis_variable.set(False)



class snapshot_GUI:
  def __init__(self, app, use_analysis):
    self.app = app
    self.parameters = plot_parameters()
    self.parameters.use_analysis = use_analysis
    app.title("Snapshot GUI")
    app.geometry('450x200+10+10')
    
    plt.ion()
    plt.close('all')
    
    self.parameters.cs = 1
    
    # find sdf files and count
    self.parameters.pathname = os.path.abspath(os.getcwd())
    runs = glob.glob1(self.parameters.pathname,"*.sdf")
    RunCounter = len(runs)
    run_array = np.zeros(RunCounter)
    for ir in range(0, RunCounter):
      run_name = runs[ir]
      run_num = int(run_name[:-4])
      run_array[ir] = run_num
    run_array = sorted(run_array)
    self.parameters.istart = int(run_array[0])
    self.parameters.iend = int(run_array[-1])
    self.parameters.sdf_num = self.parameters.istart

    # initial data import, needed for variable selection combo box
    dat = isdf.use_sdf(self.parameters.sdf_num, self.parameters.pathname,
        use_analysis = self.parameters.use_analysis,
        istart = self.parameters.istart)
    
    # create empty figures
    self.fig = plt.figure(num=1, figsize=(7.3,6), facecolor='white')
    move_figure(self.fig, 700, 10)
    self.ax1 = plt.axes()
    #self.ax1.set_aspect('equal', 'box')
    setattr(self.ax1, 'cbar', 'None')

    self.fig2 = plt.figure(num=2, figsize=(6,6), facecolor='white')
    move_figure(self.fig2, 10, 1000)
    self.ax2 = plt.axes()
    self.ax3 = op.empty_lineout(self.fig2, self.ax2)
    
    # slider - time
    self.label_slider1 = tk.Label(app, text = "Select sdf number:")
    self.label_slider1.grid(column=0, row=0)

    self.slider1 = tk.Scale(app, from_=self.parameters.istart, to=self.parameters.iend, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.callbackFunc,
                            length  = 300, resolution = 1.0)
    self.slider1.grid(column=1, row=0)
    self.slider1.set(self.parameters.istart)

    # Combo box - variable
    self.labelTop_combo1 = tk.Label(app, text="Select variable:")
    self.labelTop_combo1.grid(column=0, row=1)

    self.combo1 = ttk.Combobox(app, values=dat.variables)
    self.combo1.grid(column=1, row=1)
    self.combo1.current(0)
  
    # check box - grid
    self.grid_variable = tk.BooleanVar(app)
    self.grid_button = tk.Checkbutton(app, text="grid", variable=self.grid_variable,
                                      onvalue=True, offvalue=False)
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
                                       variable=self.anisotropies_variable,
                                       onvalue=True, offvalue=False)
    self.anisotropies_button.deselect()
    self.anisotropies_button.grid(column=0, row=4)
    
    # check box - Logarithm
    self.log_variable = tk.BooleanVar(app)
    self.log_button = tk.Checkbutton(app, text="Take log_10",
                                       variable=self.log_variable,
                                       onvalue=True, offvalue=False)
    self.log_button.deselect()
    self.log_button.grid(column=0, row=5)
  
    # button - save fig as pdf
    self.print_button = tk.Button(app, text="Save .pdf", command=self.save_pdf)
    self.print_button.grid(column=1, row=2)
    
    # button - reset button
    self.reset_button = tk.Button(app, text="Reset zoom")
    self.reset_button.grid(column=1, row=3)
    self.reset_axis_variable = tk.BooleanVar(app)
    self.reset_axis_variable.set(True)
    
    # button - exit
    self.exit_button = tk.Button(app, text="Exit", command=self.exit_gui)
    self.exit_button.grid(column=1, row=4)
    
    # button - save video
    self.video_button = tk.Button(app, text="Save video", command=self.save_video)
    self.video_button.grid(column=1, row=5)
    
    self.app.bind('<Left>', self.leftKey)
    self.app.bind('<Right>', self.rightKey)
    self.combo1.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.grid_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.polar_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.anisotropies_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.log_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.reset_button.bind("<Button-1>", self.callbackFunc1)
    
  def callbackFunc1(self, event):
    self.reset_axis_variable.set(True)
    self.callbackFunc(event)

  def callbackFunc(self, event):
    self.parameters.grid_boolean = self.grid_variable.get()
    self.parameters.use_polar = self.polar_variable.get()
    self.parameters.sdf_num = self.slider1.get()
    self.parameters.var_name = self.combo1.get()
    self.parameters.reset_axis = self.reset_axis_variable.get()
    self.parameters.view_anisotropies = self.anisotropies_variable.get()
    self.parameters.use_log = self.log_variable.get()
    
    op.data_and_plot(self.parameters.sdf_num, self.fig, self.ax1, self.fig2, self.ax2, self.ax3, self.parameters)

    self.reset_axis_variable.set(False)

  def save_pdf(self):
    self.parameters.var_name = self.combo1.get()
    pdf_name = self.parameters.var_name + '_' + 'SDF_ {0:04d}'.format(self.parameters.sdf_num) + '.pdf'
    self.fig.savefig(pdf_name)
  
  def save_video(self):
    self.parameters.reset_axis = self.reset_axis_variable.get()
    
    filename1 = self.parameters.var_name + '.mp4'
    animation = ani.FuncAnimation(self.fig, op.data_and_plot, frames=range(self.parameters.istart, self.parameters.iend+1), fargs=(self.fig, self.ax1, self.fig2, self.ax2, self.ax3, self.parameters), repeat=False)

    writer = ani.FFMpegWriter(fps=24, bitrate=2e6)
    animation.save(filename1, writer=writer)
  
  def exit_gui(self):
    plt.close(self.fig)
    plt.close(self.fig2)
    self.app.destroy()
    sys.exit('GUI exit from button press')
  
  def leftKey(self, event):
    sdf_num = self.slider1.get()
    self.slider1.set(sdf_num - 1)

  def rightKey(self, event):
    sdf_num = self.slider1.get()
    self.slider1.set(sdf_num + 1)



class plot_parameters:
  def __init__(self):
    self.sdf_num = 0
    self.use_analysis = False
    self.cs = 0
    self.pathname = 'None'
    self.istart = 0
    self.iend = 0
    self.grid_boolean = False
    self.use_polar = False
    self.var_name = 'None'
    self.reset_axis = False
    self.view_anisotropies = False
    self.use_log = False



# Make a figure class



def main():
	"""
	"""



if __name__ == "__main__":
	main()
