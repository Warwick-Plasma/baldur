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


global qe_si
qe_si = 1.6021766208e-19

global c_si
c_si = 299792458.0

global me_si
me_si = 9.10938291e-31


def options():
    data_struct = op.data_structure()
    parameters = op.plot_parameters()
    parameters.use_analysis = True

    # find sdf files and count
    parameters.pathname = os.path.abspath(os.getcwd())
    runs = glob.glob1(parameters.pathname,"*.sdf")
    parameters.sdf_num = 114

    parameters.cross_section = 269
    # initial data import, needed for variable selection combo box
    data = isdf.use_sdf(parameters.sdf_num, parameters.pathname, use_analysis = parameters.use_analysis, istart = parameters.istart)

    num_bursts = len(data.bursts)
    particles_per_path = data.Particles_per_Path

    name = "Burst1"
    name_var = "Energy_Remaining"
    name_grid = "Cell_i"
    npart = data.Particles_per_Path.data[0,:]

    beam = getattr(data, name)
    beam_energy = getattr(data, name + '_' + name_var)
    beam_icell = getattr(data, name + '_' + name_grid)
    nrays = len(beam.data)

    saved_energies = []
    saved_nepart = []
    for iray in range(0, nrays):
      print_string = 'Processing ray {:4d}'.format(iray+1) \
                   + ' of {:4d}'.format(nrays) + '   '
      sys.stdout.write('\r' + print_string)
      sys.stdout.flush()

      icell_ray = beam_icell.data[iray].data
      ii = np.asarray((np.where(icell_ray == parameters.cross_section)))[0]
      c_ray = beam_energy.data[iray].data / npart[iray] / qe_si / 1000.0 #convert to keV
      npart_ray = npart[iray]
      for i in ii:
        if i:
          saved_energies.append(c_ray[i])
          saved_nepart.append(npart_ray)

    num_bins = 100
    plt.hist(saved_energies, num_bins, weights=saved_nepart)
    #plt.plot(saved_energies,saved_nepart)
    plt.xlabel('Electron Kinetic Energy, (keV)')
    plt.ylabel('Particle Number')

    temperature_kev = 30.0
    cut_off_kev = 300.0
    kT = 1.0e3 * qe_si * temperature_kev
    KE2 = np.linspace(0.0, cut_off_kev * qe_si * 1000.0, num=num_bins)
    prob = at.pdf_maxwellian(kT, KE2)
    plt.plot(KE2 / qe_si / 1000.0, prob / np.sum(prob) * np.sum(saved_nepart), "--", label="30keV")

    temperature_kev = 40.0
    cut_off_kev = 300.0
    kT = 1.0e3 * qe_si * temperature_kev
    KE2 = np.linspace(0.0, cut_off_kev * qe_si * 1000.0, num=num_bins)
    prob = at.pdf_maxwellian(kT, KE2)
    plt.plot(KE2 / qe_si / 1000.0, prob / np.sum(prob) * np.sum(saved_nepart), "--", label="40keV")

    temperature_kev = 50.0
    cut_off_kev = 300.0
    kT = 1.0e3 * qe_si * temperature_kev
    KE2 = np.linspace(0.0, cut_off_kev * qe_si * 1000.0, num=num_bins)
    prob = at.pdf_maxwellian(kT, KE2)
    plt.plot(KE2 / qe_si / 1000.0, prob / np.sum(prob) * np.sum(saved_nepart), "--", label="50keV")

    temperature_kev = 60.0
    cut_off_kev = 300.0
    kT = 1.0e3 * qe_si * temperature_kev
    KE2 = np.linspace(0.0, cut_off_kev * qe_si * 1000.0, num=num_bins)
    prob = at.pdf_maxwellian(kT, KE2)
    plt.plot(KE2 / qe_si / 1000.0, prob / np.sum(prob) * np.sum(saved_nepart), "--", label="60keV")

    temperature_kev = 80.0
    cut_off_kev = 300.0
    kT = 1.0e3 * qe_si * temperature_kev
    KE2 = np.linspace(0.0, cut_off_kev * qe_si * 1000.0, num=num_bins)
    prob = at.pdf_maxwellian(kT, KE2)
    plt.plot(KE2 / qe_si / 1000.0, prob / np.sum(prob) * np.sum(saved_nepart), "--", label="80keV")

    plt.legend(loc="upper right")

    plt.show()



def main(argv):
  """
  """
  options()



if __name__ == "__main__":
  main(sys.argv)
