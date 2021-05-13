# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.animation as ani
import numpy as np
import glob
import sys
import os
import tkinter as tk
from tkinter import ttk
import odin_plot as op
import import_sdf as isdf
plt.switch_backend('TkAgg')



def options(*args, **kwargs):
  """This is the primary function to call when using baldur. It can be run
  with several options defined by the kwargs. It is a wrapper for other
  functions but does create the Tkinter object tk.Tk() which allows for
  passing of button information.
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


def sdf_counter(runs, user_istart, user_iend):
    """This function counts the total number of sdf files in the current
    directory to bound the slider.
    """
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
  """This class creates plots which require data from all of the sdf files.

  [__init__] is called on creation of the class and defines when the menu's
  controls will call the other functions.
  [callbackFunc1] resets the grid and calls the plot update function
  [callbackFunc] updates the plot with the new slider values, variable choice
  etc.
  """
  def __init__(self, app, use_analysis, user_istart, user_iend):
    """This function is only called on creation of the class and initialises
    the tkinter menu. First setup is done including importing data needed for
    initialising menu controls from sdf files, then the figures are
    initialised to be populated later, the tkinter objects are setup.
    """
    self.app = app
    self.parameters = op.plot_parameters()
    self.data_struct = op.data_structure()
    self.parameters.use_analysis = use_analysis
    app.title("Time history GUI")
    app.geometry('450x200+10+10')

    plt.ion()
    plt.close('all')

    self.cross_section = 1

    # find sdf files and count
    self.parameters.pathname = os.path.abspath(os.getcwd())
    runs = glob.glob1(self.parameters.pathname,"*.sdf")
    self.parameters.istart, self.parameters.iend, self.parameters.sdf_num \
        = sdf_counter(runs, user_istart, user_iend)

    # initial data import, needed for variable selection combo box
    data = isdf.use_sdf(self.parameters.sdf_num, self.parameters.pathname,
        use_analysis = self.parameters.use_analysis,
        istart = self.parameters.istart)
    self.data_struct.data = isdf.get_data_all(data, self.parameters.istart,
        self.parameters.iend, self.parameters.pathname,
        self.parameters.use_analysis, self.cross_section)

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

    control_row = 0

    # Combo box - variable
    self.labelTop_combo1 = tk.Label(app, text = "Select variable:")
    self.labelTop_combo1.grid(column=0, row=control_row)

    self.combo1 = ttk.Combobox(app, values = data.variables)
    self.combo1.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.combo1.current(0)
    self.combo1.grid(column=1, row=control_row)
    control_row += 1

    # Combo box - time variable
    self.labelTop_combo2 = tk.Label(app, text = "Select time variable:")
    self.labelTop_combo2.grid(column=0, row=control_row)

    self.combo2 = ttk.Combobox(app, values = data.variables_time)
    self.combo2.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.combo2.current(0)
    self.combo2.grid(column=1, row=control_row)
    control_row += 1

    # Combo box - grid choice
    grid_list = ['default', 'initial', 'cell number']
    self.labelTop_combo3 = tk.Label(app, text = "Select a grid:")
    self.labelTop_combo3.grid(column=0, row=control_row)

    self.combo3 = ttk.Combobox(app, values = grid_list)
    self.combo3.bind("<<ComboboxSelected>>", self.callbackFunc1)
    self.combo3.current(0)
    self.combo3.grid(column=1, row=control_row)
    control_row += 1

    # slider - scale colour
    self.label_slider1 = tk.Label(app, text = "Colourbar scale:")
    self.label_slider1.grid(column=0, row=control_row)

    self.slider1 = tk.Scale(app, from_ = -10, to = 0, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.callbackFunc,
                            length  = 200, resolution = 0.01)
    self.slider1.set(0)
    self.slider1.grid(column=1, row=control_row)
    control_row += 1

    # button - reset button
    self.reset_button = tk.Button(app, text="Reset zoom")
    self.reset_axis_variable = tk.BooleanVar(app)
    self.reset_axis_variable.set(True)
    self.reset_button.grid(column=0, row=control_row)
    control_row += 1

    # button - exit
    self.exit_button = tk.Button(app, text="Exit", command=self.exit_gui)
    self.exit_button.grid(column=0, row=control_row)
    control_row += 1

    self.reset_button.bind("<Button-1>", self.callbackFunc1)

  def callbackFunc1(self, event):
    """Reset grid and update plot
    """
    self.reset_axis_variable.set(True)
    self.callbackFunc(event)

  def callbackFunc(self, event):
    """Update 1D and 2d plots with values given by tkinter controls
    """
    self.parameters.var_name = self.combo1.get()
    self.parameters.var_name2 = self.combo2.get()
    self.parameters.grid_choice = self.combo3.get()
    self.parameters.cbar_colour_scale = self.slider1.get()
    self.parameters.reset_axis = self.reset_axis_variable.get()

    op.time_history(self.fig, self.ax1, self.cax1,
                    parameters = self.parameters, data = self.data_struct)
    op.time_history_lineout(self.fig2, self.ax2, self.ax3,
                            parameters = self.parameters, data = self.data_struct)

    self.reset_axis_variable.set(False)

  def exit_gui(self):
    """Closes all created figure windows and stops code
    """
    plt.close(self.fig)
    plt.close(self.fig2)
    self.app.destroy()
    sys.exit('GUI exit from button press')



class snapshot_GUI:
  """This class creates plots which require data from a single sdf file.

  [__init__] is called on creation of the class and defines when the menu's
  controls will call the other functions.
  [callbackFunc1] resets the grid and calls the plot update function
  [callbackFunc] updates the plot with the new slider values, variable
  choice etc.
  [left_key] update slider with arrow key
  [right_key] update slider with arrow key
  [up_key] move slider to final sdf file
  [down_key] move slider to first sdf file

  title is self explanatory:
  [save_pdf]
  [save_video]
  [exit_gui]
  """
  def __init__(self, app, use_analysis, user_istart, user_iend):
    self.app = app
    self.parameters = op.plot_parameters()
    self.parameters.use_analysis = use_analysis
    app.title("Snapshot GUI")
    app.geometry('500x400+10+10')

    plt.ion()
    plt.close('all')
    user_istart = False
    user_iend = False

    # find sdf files and count
    seperator = "/"
    self.parameters.pathname = os.path.abspath(os.getcwd())
    self.parameters.parent_dir = self.parameters.pathname.split(seperator)[:-1]
    self.parameters.parent_dir = seperator.join(self.parameters.parent_dir)
    self.parameters.dir_list[0] = self.parameters.pathname.split(seperator)[-1]
    self.parameters.dir_list[1] = self.parameters.pathname.split(seperator)[-1]
    runs = glob.glob1(self.parameters.pathname,"*.sdf")
    if len(runs) < 1:
      self.parameters.dir_list = glob.glob1(self.parameters.pathname,"*")
      self.parameters.dir_list.sort()
      self.parameters.num_dir = len(self.parameters.dir_list)
      self.parameters.parent_dir = self.parameters.pathname
      self.parameters.pathname = os.path.abspath(os.getcwd()) + '/' \
          + self.parameters.dir_list[0]
      runs = glob.glob1(self.parameters.pathname,"*.sdf")
      if len(runs) < 1:
        self.app.destroy()
        sys.exit('Failed to find .sdf files')
      if self.parameters.num_dir == 1:
        self.parameters.num_dir = 2
        self.parameters.dir_list = self.parameters.dir_list \
            + self.parameters.dir_list

    self.parameters.istart, self.parameters.iend, self.parameters.sdf_num = \
        sdf_counter(runs, user_istart, user_iend)

    # initial data import, needed for variable selection combo box
    data = isdf.use_sdf(self.parameters.sdf_num, self.parameters.pathname,
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

    control_row = 0

    # slider - time
    self.label_slider1 = tk.Label(app, text = "Select sdf number:")
    self.label_slider1.grid(column=0, row=control_row)

    self.slider1 = tk.Scale(app, from_=self.parameters.istart,
                            to=self.parameters.iend, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.callbackFunc,
                            length = 300, resolution = 1.0)
    self.slider1.grid(column=1, row=control_row)
    self.slider1.set(self.parameters.istart)
    control_row += 1

    # slider - time for comparison data
    self.label_slider2 = tk.Label(app, text = "Offset comparison file:")
    self.label_slider2.grid(column=0, row=control_row)
    self.label_slider2.grid_remove()

    self.slider2 = tk.Scale(app, from_=self.parameters.istart,
                            to=self.parameters.iend, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.callbackFunc,
                            length = 300, resolution = 1.0)
    self.slider2.set(self.parameters.istart)
    self.slider2.grid(column=1, row=control_row)
    self.slider2.grid_remove()
    control_row += 1

    # Combo box - variable 1
    self.labelTop_combo1 = tk.Label(app, text="Select variable 1:")
    self.labelTop_combo1.grid(column=0, row=control_row)

    self.combo1 = ttk.Combobox(app, values=data.variables)
    self.combo1.grid(column=1, row=control_row)
    self.combo1.current(0)
    control_row += 1

    # Combo box - variable 2
    self.labelTop_combo2 = tk.Label(app, text="Select variable 2:")
    self.labelTop_combo2.grid(column=0, row=control_row)
    self.labelTop_combo2.grid_remove()

    self.combo2 = ttk.Combobox(app, values=data.variables)
    self.combo2.grid(column=1, row=control_row)
    self.combo2.current(0)
    self.combo2.grid_remove()
    control_row += 1

    # check box - grid
    self.grid_variable = tk.BooleanVar(app)
    self.grid_button = tk.Checkbutton(app, text="grid",
                                      variable=self.grid_variable,
                                      onvalue=True, offvalue=False)
    self.grid_button.deselect()
    self.grid_button.grid(column=0, row=control_row)

    # button - save fig as pdf
    self.print_button = tk.Button(app, text="Save .pdf", command=self.save_pdf)
    self.print_button.grid(column=1, row=control_row)
    control_row += 1

    # check box - polar coordinates
    self.polar_variable = tk.BooleanVar(app)
    self.polar_button = tk.Checkbutton(app, text="polar coordinates",
                                       variable = self.polar_variable,
                                       onvalue=True, offvalue=False)
    self.polar_button.deselect()
    self.polar_button.grid(column=0, row=control_row)

    # button - reset button
    self.reset_button = tk.Button(app, text="Reset zoom")
    self.reset_axis_variable = tk.BooleanVar(app)
    self.reset_axis_variable.set(True)
    self.reset_button.grid(column=1, row=control_row)
    control_row += 1

    # check box - anisotropies
    self.anisotropies_variable = tk.BooleanVar(app)
    self.anisotropies_button = \
        tk.Checkbutton(app, text="View anisotropies",
                       variable=self.anisotropies_variable,
                       onvalue=True, offvalue=False)
    self.anisotropies_button.deselect()
    self.anisotropies_button.grid(column=0, row=control_row)

    # button - exit
    self.exit_button = tk.Button(app, text="Exit", command=self.exit_gui)
    self.exit_button.grid(column=1, row=control_row)
    control_row += 1

    # check box - Logarithm
    self.log_variable = tk.BooleanVar(app)
    self.log_button = tk.Checkbutton(app, text="Take log_10",
                                     variable=self.log_variable,
                                     onvalue=True, offvalue=False)
    self.log_button.deselect()
    self.log_button.grid(column=0, row=control_row)

    # button - save video
    self.video_button = tk.Button(app, text="Save video",
                                  command=self.save_video)
    self.video_button.grid(column=1, row=control_row)
    control_row += 1

    # Combo box - surface tracking
    self.label_combo_surf = tk.Label(app, text="Select a surface to track:")
    self.label_combo_surf.grid(column=0, row=control_row)

    data.track_surfaces.insert(0,'None')
    self.combo_surf = ttk.Combobox(app, values=data.track_surfaces)
    self.combo_surf.current(0)
    self.combo_surf.grid(column=1, row=control_row)
    control_row += 1

    # Entry - Change cbar scale max
    self.apply_scale_max = tk.BooleanVar(app)
    self.scale_max_check = tk.Checkbutton(app, text="Apply max scaling",
                                          variable=self.apply_scale_max,
                                          onvalue=True, offvalue=False)
    self.scale_max_check.deselect()
    self.scale_max_check.grid(column=0, row=control_row)

    self.entry_scale_max = tk.Entry(app)
    self.entry_scale_max.insert(0, "1.0")
    self.entry_scale_max.grid(column=1, row=control_row)
    control_row += 1

    # Entry - Change cbar scale min
    self.apply_scale_min = tk.BooleanVar(app)
    self.scale_min_check = tk.Checkbutton(app, text="Apply min scaling",
                                          variable=self.apply_scale_min,
                                          onvalue=True, offvalue=False)
    self.scale_min_check.deselect()
    self.scale_min_check.grid(column=0, row=control_row)

    self.entry_scale_min = tk.Entry(app)
    self.entry_scale_min.insert(0, "0.0")
    self.entry_scale_min.grid(column=1, row=control_row)
    control_row += 1

    # check box - plot rays?
    self.rays_variable = tk.BooleanVar(app)
    self.rays_button = tk.Checkbutton(app, text="Show rays",
                                      variable=self.rays_variable,
                                      onvalue=True, offvalue=False)
    self.rays_button.deselect()
    self.rays_button.grid(column=0, row=control_row)

    # check box - plot all rays?
    self.all_rays_variable = tk.BooleanVar(app)
    self.all_rays_button = tk.Checkbutton(app, text="All rays?",
                                      variable=self.all_rays_variable,
                                      onvalue=True, offvalue=False)
    self.all_rays_button.deselect()
    self.all_rays_button.grid(column=1, row=control_row)
    self.all_rays_button.grid_remove()
    control_row += 1

    # slider - ray selecter
    self.label_slider3 = tk.Label(app, text = "Which rays to plot:")
    self.label_slider3.grid(column=0, row=control_row)
    self.label_slider3.grid_remove()

    self.slider3 = tk.Scale(app, from_=1,
                            to=1000, tickinterval=100,
                            orient=tk.HORIZONTAL, command=self.update_ray,
                            length = 300, resolution = 1.0)
    self.slider3.set(1)
    self.slider3.grid(column=1, row=control_row)
    self.slider3.grid_remove()
    control_row += 1

    # Entry - Plot second file
    self.parameters.apply_comparison = [False] * self.parameters.num_dir
    self.parameters.apply_comparison[0] = True
    self.parameters.entry_comparison = [None] * self.parameters.num_dir
    self.parameters.entry_comparison[0] = self.parameters.pathname
    self.apply_comparison = [None] * (self.parameters.num_dir - 1)
    self.comparison_check = [None] * (self.parameters.num_dir - 1)
    self.entry_comparison = [None] * (self.parameters.num_dir - 1)
    for num in range(0, self.parameters.num_dir-1):
      self.apply_comparison[num] = tk.BooleanVar(app)
      self.comparison_check[num] = tk.Checkbutton(app, text="Apply Comparison" +
                                                  " from file ->",
                                                  variable=self.apply_comparison[num],
                                                  onvalue=True, offvalue=False)
      self.comparison_check[num].deselect()
      self.comparison_check[num].grid(column=0, row=control_row)

      self.entry_comparison[num] = tk.Entry(app)
      self.entry_comparison[num].insert(0, self.parameters.parent_dir + '/'
          + self.parameters.dir_list[num+1])
      self.entry_comparison[num].grid(column=1, row=control_row)
      control_row += 1

    # Entry - Cross section
    self.y_dir_cross_section = tk.BooleanVar(app)
    self.cross_section_check = tk.Checkbutton(app, text="Cross section X/Y, location:",
                                           variable=self.y_dir_cross_section,
                                           onvalue=True, offvalue=False)
    self.cross_section_check.deselect()
    self.cross_section_check.grid(column=0, row=control_row)

    self.entry_cross_section = tk.Entry(app)
    self.entry_cross_section.insert(0, "1")
    self.entry_cross_section.grid(column=1, row=control_row)
    control_row += 1

    # check box - show legend
    self.legend_variable = tk.BooleanVar(app)
    self.legend_button = tk.Checkbutton(app, text="Show legend",
                                      variable=self.legend_variable,
                                      onvalue=True, offvalue=False)
    self.legend_button.deselect()
    self.legend_button.grid(column=0, row=control_row)

    self.parameters.line_labels = [None] * self.parameters.num_dir
    self.entry_line = [None] * self.parameters.num_dir
    for num in range(0, self.parameters.num_dir):
      self.parameters.line_labels[num] = self.parameters.dir_list[num]
      self.entry_line[num] = tk.Entry(app)
      self.entry_line[num].insert(0, self.parameters.dir_list[num])
      self.entry_line[num].grid(column=1, row=control_row)
      control_row += 1

    # Bindings
    self.app.bind('<Left>', self.left_key)
    self.app.bind('<Right>', self.right_key)
    self.app.bind('<Up>', self.up_key)
    self.app.bind('<Down>', self.down_key)
    self.combo1.bind("<<ComboboxSelected>>", self.show_variable_combo)
    self.combo2.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.grid_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.polar_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.anisotropies_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.log_button.bind("<ButtonRelease-1>", self.callbackFunc1)
    self.reset_button.bind("<Button-1>", self.callbackFunc1)
    self.combo_surf.bind("<<ComboboxSelected>>", self.callbackFunc)
    self.scale_max_check.bind("<ButtonRelease-1>", self.callbackFunc)
    self.rays_button.bind("<ButtonRelease-1>", self.show_ray_slider)
    self.all_rays_button.bind("<ButtonRelease-1>", self.show_ray_slider)
    for num in range(0, self.parameters.num_dir - 1):
      self.comparison_check[num].bind("<ButtonRelease-1>", self.hide_slider)
    self.legend_button.bind("<ButtonRelease-1>", self.callbackFunc)
    self.cross_section_check.bind("<ButtonRelease-1>", self.callbackFunc1)

  def callbackFunc(self, event):
    """Update 1D and 2d plots with values given by tkinter controls which are
    saved in class 'parameters'
    """
    self.parameters.first_call = False
    self.parameters.grid_boolean = self.grid_variable.get()
    self.parameters.use_polar = self.polar_variable.get()
    self.parameters.sdf_num = self.slider1.get()
    self.parameters.var_name[0] = self.combo1.get()
    self.parameters.var_name[1] = self.combo2.get()
    self.parameters.reset_axis = self.reset_axis_variable.get()
    self.parameters.view_anisotropies = self.anisotropies_variable.get()
    self.parameters.use_log = self.log_variable.get()
    self.parameters.surface_name = self.combo_surf.get()
    self.parameters.apply_scale_max = self.apply_scale_max.get()
    self.parameters.scale_max = float(self.entry_scale_max.get())
    self.parameters.apply_scale_min = self.apply_scale_min.get()
    self.parameters.scale_min = float(self.entry_scale_min.get())
    self.parameters.cross_section = int(self.entry_cross_section.get())
    self.parameters.show_legend = self.legend_variable.get()
    self.parameters.line_labels[0] = self.entry_line[0].get()
    for num in range(0, self.parameters.num_dir - 1):
      self.parameters.apply_comparison[num+1] = self.apply_comparison[num].get()
      self.parameters.entry_comparison[num+1] = self.entry_comparison[num].get()
      self.parameters.line_labels[num+1] = self.entry_line[num+1].get()
    self.parameters.sdf_num2 = self.slider2.get()
    self.parameters.y_dir_cross_section = self.y_dir_cross_section.get()

    self.parameters, self.data_struct \
                    = op.data_and_plot(self.parameters.sdf_num, self.fig,
                                       self.ax1, self.cax1, self.fig2,
                                       self.ax2, self.ax3, self.parameters)

    self.reset_axis_variable.set(False)

  def show_variable_combo(self, event):
    self.parameters.var_name[0] = self.combo1.get()
    if self.parameters.var_name[0] == 'None':
      self.combo2.set('None')
      self.labelTop_combo2.grid_remove()
      self.combo2.grid_remove()
    else:
      self.labelTop_combo2.grid()
      self.combo2.grid()
    self.callbackFunc(event)

  def show_ray_slider(self, event):
    self.parameters.select_ray = self.slider3.get()
    self.parameters.plot_rays_on = self.rays_variable.get()
    self.parameters.plot_all_rays = self.all_rays_variable.get()
    if self.parameters.plot_rays_on:
      self.label_slider3.grid()
      self.slider3.grid()
      self.all_rays_button.grid()
      if (self.parameters.first_call):
        self.callbackFunc(event)
      else:
        op.wrapper_plot_light_rays(self.data_struct.data[0], self.parameters,
            self.fig, self.ax1)
        op.wrapper_plot_electron_rays(self.data_struct.data[0], self.parameters,
            self.fig, self.ax1)
    else:
      self.label_slider3.grid_remove()
      self.slider3.grid_remove()
      self.all_rays_button.grid_remove()

  def update_ray(self, event):
    self.parameters.select_ray = self.slider3.get()
    self.parameters.plot_rays_on = self.rays_variable.get()
    self.parameters.plot_all_rays = self.all_rays_variable.get()
    op.wrapper_plot_electron_rays(self.data_struct.data[0], self.parameters, self.fig,
        self.ax1)

  def hide_slider(self, event):
    show_slider_boolean = self.apply_comparison[0].get()
    if show_slider_boolean:
      print()
      print("Use arrow keys to move sliders together")

      self.slider2.set(self.slider1.get())
      self.label_slider2.grid()
      self.slider2.grid()
    else:
      self.label_slider2.grid_remove()
      self.slider2.grid_remove()
    self.callbackFunc(event)

  def callbackFunc1(self, event):
    """This function resets the grid before the plot updateing function
    """
    self.reset_axis_variable.set(True)
    self.callbackFunc(event)

  def save_pdf(self):
    """Save pdf of 2d file as default but easily changed by using fig1
    instead of fig
    """
    self.parameters.var_name = self.combo1.get()
    pdf_name = self.parameters.var_name + '_' + \
        'SDF_{0:04d}'.format(self.parameters.sdf_num) + '.pdf'
    self.fig.savefig(pdf_name)

  def save_video(self):
    """Saves a video of 2d plot as default but can be changed by changing fig
    """
    self.parameters.reset_axis = self.reset_axis_variable.get()
    self.parameters.movie = True

    filename1 = self.parameters.var_name[0] + '.gif'
    animation = ani.FuncAnimation(self.fig, op.data_and_plot,
        frames=range(self.parameters.sdf_num, self.parameters.iend+1),
        fargs=(self.fig, self.ax1, self.cax1, self.fig2, self.ax2, self.ax3,
        self.parameters))

    writer = PillowWriter(fps=24)
    animation.save(filename1, writer=writer)
    self.parameters.movie = False

  def exit_gui(self):
    """Closes all created figure windows and stops code
    """
    plt.close(self.fig)
    plt.close(self.fig2)
    self.app.destroy()
    sys.exit('GUI exit from button press')

  def left_key(self, event):
    """hotkey updates slider
    """
    sdf_num = self.slider1.get()
    sdf_num2 = self.slider2.get()
    self.slider1.set(sdf_num - 1)
    self.slider2.set(sdf_num2 - 1)

  def right_key(self, event):
    """hotkey updates slider
    """
    sdf_num = self.slider1.get()
    sdf_num2 = self.slider2.get()
    self.slider1.set(sdf_num + 1)
    self.slider2.set(sdf_num2 + 1)

  def up_key(self, event):
    """hotkey updates slider
    """
    self.slider1.set(self.parameters.iend)
    self.slider2.set(self.parameters.iend)

  def down_key(self, event):
    """hotkey updates slider
    """
    self.slider1.set(self.parameters.istart)
    self.slider2.set(self.parameters.istart)



def main(argv):
  """
  """
  if (len(argv) != 3):
    print("Provide options at command line i.e.")
    print("python menu.py 1 0")
    print("meaning use_analysis = True and time_history = False")
    sys.exit('Exit due to insufficient arguments')
  if (int(argv[1]) == 1):
    analysis = True
  else:
    analysis = False
  if (int(argv[2]) == 1):
    history = True
  else:
    history = False
  options(use_analysis = analysis, time_history = history)



if __name__ == "__main__":
  main(sys.argv)
