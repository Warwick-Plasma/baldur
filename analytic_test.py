# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons, TextBox
from scipy.optimize import fsolve
import import_sdf as isdf
import csv

global qe_si
qe_si = 1.6021766208e-19

global c_si
c_si = 299792458.0

global me_si
me_si = 9.10938291e-31

global rest_mass_energy
rest_mass_energy = me_si * c_si**2



def compare_penetration_function():
  """
  """
  fig = plt.figure()
  ax = plt.axes()

  csv_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_unionized_50_keV/Berger_seltzer_maeda_penetration_50kev.csv"

  with open(csv_name) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    depths = []
    pen_funcs = []
    for row in readCSV:
      depth = float(row[0])
      depths.append(depth)
      pen_func = float(row[1])
      pen_funcs.append(pen_func)

  ax.plot(depths, pen_funcs, "k-", label="Grun et al. 1957 experimental")

  path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_berger_model/0000.sdf"
  plot_penetration_function(path_name, ax, legend_label="Berger, Seltzer and Maeda model")

  path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_fox_non_rel_model/0000.sdf"
  plot_penetration_function(path_name, ax, legend_label="Fox ionized, non-rel model")

  path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_davies_non_rel_model/0000.sdf"
  plot_penetration_function(path_name, ax, legend_label="Davies unionized, non-rel model")

  #path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_hybrid_atom_scatter_ion_stopping/0000.sdf"
  #plot_penetration_function(path_name, ax, legend_label="Ion stopping, atom scattering hyrid model")

  #path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_hybrid_ion_scatter_atom_stopping/0000.sdf"
  #plot_penetration_function(path_name, ax, legend_label="Atom stopping, ion scattering hyrid model")

  path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_robinson_ionized_model/0000.sdf"
  plot_penetration_function(path_name, ax, legend_label="Robinson ionized model")

  path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_robinson_unionized_model/0000.sdf"
  plot_penetration_function(path_name, ax, legend_label="Robinson unionized model")

  path_name = "/home/phrjmt/Odin/Data_electron_dep_berger_seltzer_maeda_1969_charge_state_half/0000.sdf"
  plot_penetration_function(path_name, ax, legend_label="Davies/Fox q=0.5, non-rel model")

  ax.set_xlabel("Penetration Depth Normalised to Mean Range (Areal Density)")
  ax.set_ylabel("Dose Function Normalised by Electron Energy and Mean Range")
  ax.set_title("50keV Mono-energetic Hot Electron Stopping into Cold Low Pressure Nitrogen, Experimental Verification")
  ax.legend()

  plt.show()



def plot_penetration_function(path_name, ax, legend_label):
  """
  """
  dat = sh.getdata(path_name)
  hot_electron_energy_MeV = 50.0e-3
  mean_range = 4.92e-3

  grid_mid = dat.Grid_Grid_mid.data
  xc = grid_mid[0]
  yc = grid_mid[1]

  electron_dep = dat.Fluid_Energy_deposited_hot_electron.data
  nepart = np.sum(dat.Particles_per_Path.data)
  inverse_electron_charge = 6.242e+18
  electron_erg_per_vol = electron_dep / nepart * inverse_electron_charge / 1.0e6 / 1.0e6 # MeV -> J and m^-3 -> cm^-3 (mass cancels)
  units_new = "MeV/cm$^3$/source particle"

  dr = np.zeros(np.shape(xc))
  for i in range(len(dr)-1):
    dr[i,:] = xc[i+1,:] - xc[i,:]

  penetration_function = np.sum((2.0*np.pi*electron_erg_per_vol*xc*dr*1.0e4), axis=0)*1000.0
  pen_func_normalised = penetration_function * mean_range / hot_electron_energy_MeV

  dr = np.zeros(np.shape(yc[0,:]))
  for i in range(len(dr)-1):
    dr[i] = yc[0,i+1] - yc[0,i]
  rho = dat.Fluid_Rho.data[0,:]
  rhor = rho * dr
  rhor_cumulative = np.cumsum(rhor, axis=0)
  rhor_gcm = rhor_cumulative * 0.1
  rhor_gcm_normalised = rhor_gcm / mean_range

  ax.plot(rhor_gcm_normalised, pen_func_normalised[::-1], ".", label=legend_label)



def plot_maxwellian(*args, **kwargs):
  """
  """
  sample_size = kwargs.get('sample_size', 1000)
  temperature_kev = kwargs.get('temperature_kev', 50.0)
  kT = 1.0e3 * qe_si * temperature_kev
  num_bins = 100

  KE = direct_sample_maxwellian(sample_size, kT)

  KE2 = np.linspace(np.min(KE), np.max(KE), num=num_bins)
  prob = pdf_maxwellian(kT, KE2)

  KE = KE / rest_mass_energy
  KE2 = KE2 / rest_mass_energy

  weights = np.ones_like(KE) / len(KE)
  plt.hist(KE, num_bins, weights=weights)
  plt.style.use('ggplot')

  plt.plot(KE2, prob / np.sum(prob), lw=3)

  plt.xlabel('Electron Kinetic Energy, units of rest mass')
  plt.ylabel('Probability')

  plt.show()


def pdf_maxwellian(kT, KE):
  pdf = 2.0 * np.sqrt(KE / np.pi) * (kT)**(-1.5) * np.exp(-KE / kT)

  return pdf

def direct_sample_maxwellian(sample_size, kT):

  rand_num = np.random.random((3, sample_size))
  KE = kT * (-np.log(rand_num[0,:]) - np.log(rand_num[1,:])
      * np.cos(np.pi / 2.0 * rand_num[2,:])**2)

  return KE



def plot_multi():
  """
  """
  sdf_num = 10
  search_string = str(sdf_num).zfill(4)+'.sdf'
  pathname = os.path.abspath(os.getcwd())
  print('From directory:')
  print(pathname)

  sub_dir = sorted(glob.glob1(pathname, '*'))

  plt.figure()
  fig_labels = ["arithmetic", "harmonic $\delta = 0.1$",
      "harmonic $\delta = 0.01$", "harmonic $\delta = 0.001$"]
  fig_symbols = ["","x",".","+"]
  fig_colour = ["r","g","b","m"]
  fig_line = ["-","None","None","None"]


  for idir in range(0, len(sub_dir)):
    print('We have imported from directories:')
    print(sub_dir[idir])

    current_dir_ap = pathname + '/' + sub_dir[idir]
    dat = isdf.use_sdf(sdf_num, current_dir_ap)

    extent = np.shape(dat.Fluid_Temperature_electron.data)
    cs = int(extent[1] * 0.5)
    temp = dat.Fluid_Temperature_electron.data[:,cs]
    x = dat.Grid_Grid_mid.data[0][:,cs]

    plt.plot(x, temp, label = fig_labels[idir], color = fig_colour[idir],
        marker = fig_symbols[idir], linestyle = fig_line[idir])

  plt.title("Non-linear heatwave travelling into cold medium")
  plt.xlabel("Distance (arbitrary)")
  plt.ylabel("Temperature (arbitrary)")
  plt.legend()

  plt.show()



def analytic_comparison(istart):
  Nt = 100
  t0 = int(Nt * istart / 10) - 1
  cs = 18

  x_sim, temp_sim = simulation_data(istart)
  x_ana, temp_ana = non_linear_heat_wave(x_sim, use_time = True, time = 1.0)

  L1, L2, Linfty = error_analysis(temp_sim, temp_ana)

  plt.figure()
  plt.plot(x_ana[:,cs], temp_ana[:,cs], 'r', linewidth=4)
  plt.plot(x_sim, temp_sim, 'kx', markersize=2)
  plt.title("Non-linear heatwave travelling into warm medium, kershaw grid")
  plt.xlabel("Distance (arbitrary)")
  plt.ylabel("Temperature (arbitrary)")

  initial_text = "harmonic\n512 time steps\n" + "L1 Norm = {:.3e}".format(L1) + "\n" + "L2 Norm = {:.3e}".format(L2) + "\n" + "L$_\infty$ Norm = {:.3e}".format(Linfty)
  axbox = plt.axes([0.15, 0.15, 0.3, 0.2])
  text_box = TextBox(axbox, '', initial=initial_text)

  plt.show()



def simulation_data(istart):
  pathname = os.path.abspath(os.getcwd())
  dat = isdf.use_sdf(istart, pathname)

  extent = np.shape(dat.Fluid_Temperature_electron.data)
  cs = int(extent[1] * 0.5)
  temp = dat.Fluid_Temperature_electron.data#[:,cs]
  x = dat.Grid_Grid_mid.data[0]#[:,cs]

  return x, temp



def non_linear_heat_wave(x, *args, **kwargs):
  use_time = kwargs.get('use_time', False)
  time = kwargs.get('time', 0.0)

  if use_time:
    Nx = np.shape(x)[0]
    Nt = np.shape(x)[1]
    T = time
    X = x
  else:
    Nx = np.shape(x)[0]#38
    #dx = (1.0 / (Nx - 2)) / 2
    #x = np.linspace(-dx, 1.0 + dx, Nx)
    Nt = 100
    t = np.linspace(0, 1, Nt)
    X, T = np.meshgrid(x, t, indexing='ij')

  kappa_0 = 0.01
  XI = (1 / kappa_0) * ((T - X) - 0.15)
  PHI = XI # initialise PHI

  func = lambda phi, xi : phi - np.exp(xi -(4.0 * np.abs(phi) + 3.0 * phi**2 + 4.0 / 3.0 * np.abs(phi)**3 + 1.0 / 4.0 * phi**4))

  for ix in range(0, Nx):
    for it in range(0, Nt):
      phi_guess = 0.4
      args = XI[ix, it]
      PHI[ix, it] = fsolve(func, phi_guess, args=args, xtol=1e-06, maxfev=500)

  temp = 1 + PHI

  end_string = '_num, '
  temp_list = [None] * Nt

  for it in range(0, Nt):
    nt = it + 1
    if nt % 5 == 0:
      temp_list[it] = "{:.5f}_num, &".format(temp[0,it])
    else:
      temp_list[it] = "{:.5f}_num, ".format(temp[0,it])

  with open('output.txt', 'w') as f:
    for it in range(0, Nt):
      nt = it + 1
      if nt % 5 == 0:
        f.write("%s\n" % temp_list[it])
      else:
        f.write("%s" % temp_list[it])

  return X, temp



def error_analysis(temp_sim, temp_ana):
  L1 = np.mean(abs(temp_sim - temp_ana))
  L2 = np.sqrt(np.mean(abs(temp_sim - temp_ana)**2))
  Linfty = np.max(np.abs(temp_sim - temp_ana))

  print("The L1 Norm is ", L1)
  print("The L2 Norm is ", L2)
  print("The Linfty Norm is ", Linfty)
  return L1, L2, Linfty



def main():
	"""
	"""


if __name__ == "__main__":
	main()
