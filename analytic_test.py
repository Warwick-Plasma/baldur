# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons
from scipy.optimize import fsolve
import import_sdf as isdf



def analytic_comparison(istart):
  Nt = 100
  t0 = int(Nt * istart / 10) - 1
  cs = 50
  
  x_sim, temp_sim = simulation_data(istart)
  x_ana, temp_ana = non_linear_heat_wave(x_sim[:,cs])

  L1, L2, Linfty = error_analysis(temp_sim[:,cs], temp_ana[:,t0])
  
  print("The L1 Norm is ", L1)
  print("The L2 Norm is ", L2)
  print("The Linfty Norm is ", Linfty)

  plt.figure()
  plt.plot(x_ana[:,t0], temp_ana[:,t0], 'r')
  plt.plot(x_sim[:,cs], temp_sim[:,cs], 'kx')
  plt.title("Non-linear heatwave travelling")
  plt.xlabel("Distance (arbitrary)")
  plt.ylabel("Temperature (arbitrary)")
  
  plt.show()



def simulation_data(istart):
  analysis = False  
  pathname = os.path.abspath(os.getcwd())
  dat = isdf.use_sdf(istart, pathname, analysis)
  
  extent = np.shape(dat.Fluid_Temperature_electron.data)
  cs = int(extent[1] * 0.5)
  temp = dat.Fluid_Temperature_electron.data#[:,cs]
  x = dat.Grid_Grid_mid.data[0]#[:,cs]
  
  return x, temp



def non_linear_heat_wave(x):
  Nx = np.shape(x)[0]
  Nt = 100
  t = np.linspace(0, 100, Nt)
  X, T = np.meshgrid(x, t, indexing='ij')
  kappa_0 = 1.0
  
  XI = (1 / kappa_0) * (T - X)
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
  return L1, L2, Linfty



def main():
	"""
	"""


if __name__ == "__main__":
	main()
