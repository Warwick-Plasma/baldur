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
import odin_plot as op
import import_sdf as isdf
import menu as mu
import analytic_test as at
plt.switch_backend('TkAgg')



def options():
    data_struct = op.data_structure()
    user_istart = 100
    user_iend = 130
    parameters = op.plot_parameters()
    parameters.use_analysis = True

    # find sdf files and count
    parameters.pathname = os.path.abspath(os.getcwd())
    runs = glob.glob1(parameters.pathname,"*.sdf")
    parameters.istart, parameters.iend, parameters.sdf_num = mu.sdf_counter(runs, user_istart, user_iend)

    parameters.cross_section = 1
    # initial data import, needed for variable selection combo box
    data = isdf.use_sdf(parameters.sdf_num, parameters.pathname, use_analysis = parameters.use_analysis, istart = parameters.istart)
    data_struct.data = isdf.get_data_all(data, parameters.istart, parameters.iend, parameters.pathname, parameters.use_analysis, parameters.cross_section)

    fig1 = plt.figure(num=1, figsize=(6,6), facecolor='white')
    mu.move_figure(fig1, 10, 1000)
    ax1 = plt.axes()
    ax2 = op.empty_lineout(fig1, ax1)
    ax2.set_visible(False)

    at.lawson_criteria(ax1)
    op.plot_thermodynamic_path(fig1, ax1, data_struct.data, parameters)



def main(argv):
  """
  """
  options()



if __name__ == "__main__":
  main(sys.argv)
