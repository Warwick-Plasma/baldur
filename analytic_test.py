# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons, TextBox
from scipy.optimize import fsolve
import import_sdf as isdf



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
