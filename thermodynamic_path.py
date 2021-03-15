# Duncan Barlow, Odin project, Warwick University, 02/21
import sdf_helper as sh
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from matplotlib.animation import FuncAnimation
import matplotlib.animation as ani
import numpy as np
import glob
import sys
import os
import tkinter as tk
from tkinter import ttk
import csv
import odin_plot as op
import import_sdf as isdf
import menu as mu
import analytic_test as at
plt.switch_backend('TkAgg')



def options():
    data_struct = op.data_structure()
    user_istart = 111
    user_iend = 130
    parameters = op.plot_parameters()
    parameters.use_analysis = True

    parameters.pathname, runs, parameters.dir_list, \
        parameters.num_dir, parameters.parent_dir = mu.find_sdf_files()

    # find sdf files and count
    parameters.istart, parameters.iend, parameters.sdf_num = mu.sdf_counter(runs, user_istart, user_iend)

    fig1 = plt.figure(num=1, figsize=(6,6), facecolor='white')
    mu.move_figure(fig1, 10, 1000)
    ax1 = plt.axes()
    ax2 = op.empty_lineout(fig1, ax1)
    ax2.set_visible(False)

    at.lawson_criteria(ax1)
    for num in range(0, parameters.num_dir):
      pathname = parameters.parent_dir + '/' + parameters.dir_list[num]
      if os.path.isdir(pathname):
        file_name = pathname+"/thermodynamic_path.csv"
        if os.path.isfile(file_name):
          with open(file_name, "r") as file:
            reader = csv.reader(file)
            x_data = reader.__next__()
            y_data = reader.__next__()
          x_label = x_data[0]
          y_label = y_data[0]

          x_data = np.array(x_data[1:]).astype(np.float)
          y_data = np.array(y_data[1:]).astype(np.float)
        else:
          parameters.cross_section = 1
          # initial data import, needed for variable selection combo box
          data = isdf.use_sdf(parameters.sdf_num, parameters.pathname, \
              use_analysis = parameters.use_analysis, istart = parameters.istart)
          data = isdf.get_data_all(data, parameters.istart, parameters.iend, \
              parameters.pathname, parameters.use_analysis, parameters.cross_section)

          var = data.Ion_Temperature_Mean_Hotspot
          y_data = var.all_time_data * var.unit_conversion
          y_label = var.name + " (" + var.units_new + ")"

          rhor = data.Areal_Density_Mean_Hotspot
          x_data = rhor.all_time_data * rhor.unit_conversion
          x_label = rhor.name + " (" + rhor.units_new + ")"

          file_name = "thermodynamic_path"
          op.save_as_csv(file_name, x_data, y_data, xlabel=x_label, ylabel=y_label)

        line_label = parameters.dir_list[num]
        op.plot_thermodynamic_path(fig1, ax1, num+1, x_data, y_data, \
            x_label, y_label, line_label)
      else:
        print()
        print("Warning: " + pathname + " is not a directory")
    ax1.legend(loc = 'upper right')
    plt.show()



def main(argv):
  """
  """
  options()



if __name__ == "__main__":
  main(sys.argv)
