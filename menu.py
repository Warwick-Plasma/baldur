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
plt.switch_backend('TkAgg')



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
  """This function creates the Tkinter object tk.Tk()
  """
  use_analysis = kwargs.get('use_analysis', False)
  user_istart = kwargs.get('user_istart', False)
  user_iend = kwargs.get('user_iend', False)
  time_history = kwargs.get('time_history', False)
  op.check_analysis(use_analysis)
  root = tk.Tk()
  if (time_history == False):
    print('Set <time_history = True>?')
    my_gui = snapshot_GUI(root, use_analysis, user_istart, user_iend)
  else:
    my_gui = time_history_GUI(root, use_analysis, user_istart, user_iend)
  root.mainloop()



def sdf_counter(runs, user_istart, user_iend):
    RunCounter = len(runs)
    run_array = np.zeros(RunCounter)
    for ir in range(0, RunCounter):
      run_name = runs[ir]
      run_num = int(run_name[:-4])
      run_array[ir] = run_num
    run_array = sorted(run_array)
    if user_istart==False:
      istart = int(run_array[0])
    else:
      istart = user_istart
    if user_iend==False:
      iend = int(run_array[-1])
    else:
      iend = user_iend
    sdf_num = istart
    
    return istart, iend, sdf_num



# Classes are required in these section to store the objects, figures and the
# various controls defined by tkinter
class time_history_GUI:
  """This class creates plots which require data from all of the time
  snapshot outputs.
  """
  def __init__(self, app, use_analysis, user_istart, user_iend):
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
    self.istart, self.iend, sdf_num = sdf_counter(runs, user_istart, user_iend)
    
    # initial data import, needed for variable selection combo box
    dat = isdf.use_sdf(sdf_num, self.pathname, use_analysis = self.use_analysis, istart = self.istart)
    self.dat = isdf.get_data_all(dat, self.istart, self.iend, self.pathname, self.use_analysis, self.cs)
    
    # create empty figures
    aspc = 1.2
    self.fig = plt.figure(num=1, figsize=(8*aspc,8), facecolor='white')
    move_figure(self.fig, 700, 10)
    self.ax1 = self.fig.add_axes([0.15, 0.14, 0.6, 0.6*aspc])
    self.cax1 = self.fig.add_axes([0.8, 0.14, 0.05, 0.6*aspc])

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
      
    op.time_history(self.dat, self.fig, self.ax1, self.cax1, var_name = var_name, cbar_upscale = cbar_upscale, grid = grid_choice, reset_axis = reset_axis)
      
    op.time_history_lineout(self.dat, self.fig2, self.ax2, self.ax3, var_name = var_name2,  use_analysis = self.use_analysis)

    self.reset_axis_variable.set(False)



class snapshot_GUI:
  def __init__(self, app, use_analysis, user_istart, user_iend):
    self.app = app
    self.parameters = plot_parameters()
    self.parameters.use_analysis = use_analysis
    app.title("Snapshot GUI")
    app.geometry('500x400+10+10')
    
    plt.ion()
    plt.close('all')
    
    self.parameters.cs = 1
    # find sdf files and count
    self.parameters.pathname = os.path.abspath(os.getcwd())
    runs = glob.glob1(self.parameters.pathname,"*.sdf")
    self.parameters.istart, self.parameters.iend, self.parameters.sdf_num = sdf_counter(runs, user_istart, user_iend)

    # initial data import, needed for variable selection combo box
    dat = isdf.use_sdf(self.parameters.sdf_num, self.parameters.pathname,
        use_analysis = self.parameters.use_analysis,
        istart = self.parameters.istart)
    
    # create empty figures
    aspc = 1.2
    self.fig = plt.figure(num=1, figsize=(8*aspc,8), facecolor='white')
    move_figure(self.fig, 700, 10)
    self.ax1 = self.fig.add_axes([0.15, 0.14, 0.6, 0.6*aspc])
    self.cax1 = self.fig.add_axes([0.8, 0.14, 0.05, 0.6*aspc])
    #self.ax1.set_aspect('equal', 'box')

    self.fig2 = plt.figure(num=2, figsize=(6,6), facecolor='white')
    move_figure(self.fig2, 10, 1050)
    self.ax2 = plt.axes()
    self.ax3 = op.empty_lineout(self.fig2, self.ax2)
    setattr(self.ax2 , "loc_cell_track", 0)
    
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

    # Combo box - surface tracking
    self.label_combo_surf = tk.Label(app, text="Select a surface to track:")
    self.label_combo_surf.grid(column=0, row=6)

    dat.track_surfaces.insert(0,'None')
    self.combo_surf = ttk.Combobox(app, values=dat.track_surfaces)
    self.combo_surf.grid(column=1, row=6)
    self.combo_surf.current(0)
    
    # Entry - Change cbar scale max
    self.apply_scale_max = tk.BooleanVar(app)
    self.scale_max_check = tk.Checkbutton(app, text="Apply max scaling", variable=self.apply_scale_max, onvalue=True, offvalue=False)
    self.scale_max_check.deselect()
    self.scale_max_check.grid(column=0, row=7)
    
    self.entry_scale_max = tk.Entry(app)
    self.entry_scale_max.insert(0, "1.0")
    self.entry_scale_max.grid(column=1, row=7)
    
    # Entry - Change cbar scale min
    self.apply_scale_min = tk.BooleanVar(app)
    self.scale_min_check = tk.Checkbutton(app, text="Apply min scaling", variable=self.apply_scale_min, onvalue=True, offvalue=False)
    self.scale_min_check.deselect()
    self.scale_min_check.grid(column=0, row=8)
    
    self.entry_scale_min = tk.Entry(app)
    self.entry_scale_min.insert(0, "0.0")
    self.entry_scale_min.grid(column=1, row=8)
    
    # check box - plot rays?
    self.rays_variable = tk.BooleanVar(app)
    self.rays_button = tk.Checkbutton(app, text="Show rays",
                                       variable=self.rays_variable,
                                       onvalue=True, offvalue=False)
    self.rays_button.deselect()
    self.rays_button.grid(column=0, row=9)
    
    # Bindings
    self.app.bind('<Left>', self.leftKey)
    self.app.bind('<Right>', self.rightKey)
    self.combo1.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.grid_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.polar_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.anisotropies_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.log_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.reset_button.bind("<Button-1>", self.callbackFunc1)
    self.combo_surf.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.scale_max_check.bind("<ButtonRelease-1>", self.callbackFunc)
    self.rays_button.bind("<ButtonRelease-1>", self.callbackFunc)

  def callbackFunc(self, event):
    self.parameters.grid_boolean = self.grid_variable.get()
    self.parameters.use_polar = self.polar_variable.get()
    self.parameters.sdf_num = self.slider1.get()
    self.parameters.var_name = self.combo1.get()
    self.parameters.reset_axis = self.reset_axis_variable.get()
    self.parameters.view_anisotropies = self.anisotropies_variable.get()
    self.parameters.use_log = self.log_variable.get()
    self.parameters.surface_name = self.combo_surf.get()
    self.parameters.apply_scale_max = self.apply_scale_max.get()
    self.parameters.scale_max = float(self.entry_scale_max.get())
    self.parameters.apply_scale_min = self.apply_scale_min.get()
    self.parameters.scale_min = float(self.entry_scale_min.get())
    self.parameters.plot_rays_on = self.rays_variable.get()
    
    op.data_and_plot(self.parameters.sdf_num, self.fig, self.ax1, self.cax1, self.fig2, self.ax2, self.ax3, self.parameters)

    self.reset_axis_variable.set(False)
    
  def callbackFunc1(self, event):
    """This function resets the grid before calling other functions
    """
    self.reset_axis_variable.set(True)
    self.callbackFunc(event)

  def save_pdf(self):
    self.parameters.var_name = self.combo1.get()
    pdf_name = self.parameters.var_name + '_' + 'SDF_{0:04d}'.format(self.parameters.sdf_num) + '.pdf'
    self.fig.savefig(pdf_name)
  
  def save_video(self):
    self.parameters.reset_axis = self.reset_axis_variable.get()
    
    filename1 = self.parameters.var_name + '.mp4'
    animation = ani.FuncAnimation(self.fig, op.data_and_plot, frames=range(self.parameters.istart, self.parameters.iend+1), fargs=(self.fig, self.ax1, self.cax1, self.fig2, self.ax2, self.ax3, self.parameters), repeat=False)

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
    self.surface_name = 'None'
    self.apply_scale_max = False
    self.scale_max = 1.0
    self.apply_scale_min = False
    self.scale_min = 0.0
    self.plot_rays_on = False
    



def main():
	"""
	"""



if __name__ == "__main__":
	main()
