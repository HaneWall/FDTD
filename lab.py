import fdtd_1d as f
import numpy as np


from fdtd_1d.constants import c0, mu0, eps0
import matplotlib.pyplot as plt


test = f.Grid(101, 4.3e-09) # creates 201 grid cells (á 4.3e-09m)

#test[0] = f.LeftSideGridBoundary()
#test[125] = f.RightSideGridBoundary()   # note that 201 grid cells means cell 200 is the last
#test[80:121] = f.NonDispersiveMedia(name='Medium4Eps', permeability=1, permittivity=4, conductivity=0)
#test[100:171] = f.NonDispersiveMedia(name='Media2Epsleft', permeability=1, permittivity=2, conductivity=5e04)
#test[75:171] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=5e04)
test[50:101] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)
#test[20] = f.SinusoidalImpulse(name='SinusMid', amplitude=1, phase_shift=0, wavelength=430.0e-09)
#test[20] = f.EnvelopeSinus(name='EnevelopedSinusMitte', wavelength=200e-09, phase_shift=0, amplitude=1, fwhm=600e-09, peak_timestep=200)
#test[50] = f.GaussianImpulse(name='a ndererImpuls', amplitude=1, peak_timestep=60, fwhm=100.e-09)
#test[34] = f.ActivatedSinus(name='sin**2acivation', wavelength=900.0e-09, carrier_wavelength=6000.0e-9, phase_shift=0, amplitude=1)
#test[96] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=6000)
test[0] = f.LeftSideMur()
test[100] = f.RightSideMur()

test.animate_timesteps(1000)
#test.get_observed_signals()

'''
def theo_2(d, eps, timestep):
    n = np.sqrt(eps)
    k0 = test.sources[0].omega / c0
    k = n * k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(-2j * k * d)
    e_inc = np.exp(1j*test.sources[0].omega * test.dt * timestep)
    amplitude_inc = np.imag(e_inc)
    e_tr = e_inc * (2/(n+1)) * (2*n / (n+1)) * (1/(1-q)) * np.exp(-1j*(k+k0)*d)
    amplitude_tr = np.imag(e_tr)/2
    return amplitude_tr
# theoriegraphen erstellen


data_theo_1 = []
data_theo_2 = []
for ts in range(6000, 6300):
    data_theo_2.append(theo_2(31*test.dz, 4, ts))

x = np.arange(6000, 6300, 1)
y = np.arange(0, 6300, 1)
z = np.arange(0, 6300, 1)

fig, axes = plt.subplots(1,1)
#axes.plot(x, data_theo_1, label='theo_1', linestyle='dotted')
axes.plot(x, data_theo_2, label='theo_2')
axes.plot(y, test.local_observers[0].hard_save, label='hard_save')
axes.plot(z, test.local_observers[0].amplitude * np.cos(test.sources[0].omega * test.dt * y + test.local_observers[0].phase), label='reconstructed', linestyle='dashed')
#axes.plot(x, np.sin(test.sources[0].omega*test.dt * x), label='src')
axes.plot([test.local_observers[0].first_timestep, test.local_observers[0].second_timestep], test.local_observers[0].observedE, linestyle='None', marker='o')
axes.legend(loc='best')
plt.show()
'''