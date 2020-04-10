import fdtd_1d as f
import numpy as np


from fdtd_1d.constants import c0, mu0, eps0
import matplotlib.pyplot as plt

#   Spielwiese

test = f.Grid(126, 4.3e-09) # creates 201 grid cells (á 4.3e-09m)

test[0] = f.LeftSideGridBoundary()
test[125] = f.RightSideGridBoundary()   # note that 201 grid cells means cell 200 is the last
#test[80:121] = f.NonDispersiveMedia(name='Medium4Eps', permeability=1, permittivity=4, conductivity=0)
#test[100:171] = f.NonDispersiveMedia(name='Media2Epsleft', permeability=1, permittivity=2, conductivity=5e04)
#test[75:171] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=5e04)
test[50:81] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)
#test[20] = f.SinusoidalImpulse(name='SinusMid', amplitude=1, phase_shift=0, wavelength=230.0e-09)
#test[35] = f.EnvelopeSinus(name='EnevelopedSinusMitte', wavelength=100e-09, phase_shift=0, amplitude=1, fwhm=500e-09, peak_timestep=200)
#test[30] = f.GaussianImpulse(name='andererImpuls', amplitude=1, peak_timestep=60, fwhm=120.e-09)
test[20] = f.ActivatedSinus(name='sin**2acivation', wavelength=560.0e-09, carrier_wavelength= 6000.0e-9, phase_shift=0, amplitude=1)
test[110] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=1000)
test.run_timesteps(1300)
test.get_observed_signals()

def theoretical_dielectric_disc(d, eps, timestep):      # TODO somehow wrong
    n = np.sqrt(eps)
    k0 = test.sources[0].omega / c0
    k = n * k0
    q = ((n-1)**2) / ((n+1)**2) * np.exp(2j*k*d)
    e_inc = np.exp(1j*test.sources[0].omega * test.dt * timestep) * np.exp(1j*30*test.dz * k0)
    e_tr = e_inc * (2/(n+1)) * (2*n / (n+1))  * (1 / (1-q)) * np.exp(1j*k*d)
    amplitude = np.imag(e_tr)
    return amplitude

# theoriegraphen erstellen

data_theo = []
for ts in range(1000, 1200):
    data_theo.append(theoretical_dielectric_disc(31*test.dz, 4, ts))

x = np.arange(1000, 1200, 1)

fig, axes = plt.subplots(1,1)
axes.plot(x, data_theo, label='theo')
axes.plot(x, test.local_observers[0].amplitude * np.cos(test.sources[0].omega * test.dt * x + test.local_observers[0].phase), label='measured')
axes.legend(loc='best')
plt.show()