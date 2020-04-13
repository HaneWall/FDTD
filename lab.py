import fdtd_1d as f
import numpy as np
from fdtd_1d.constants import c0
import matplotlib.pyplot as plt


test = f.Grid(101, 4.3e-09) # creates 201 grid cells (á 4.3e-09m)

#test[0] = f.LeftSideGridBoundary()
#test[125] = f.RightSideGridBoundary()   # note that 201 grid cells means cell 200 is the last
#test[80:121] = f.NonDispersiveMedia(name='Medium4Eps', permeability=1, permittivity=4, conductivity=0)
#test[100:171] = f.NonDispersiveMedia(name='Media2Epsleft', permeability=1, permittivity=2, conductivity=5e04)
#test[75:171] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=5e04)
test[50:81] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)
#test[20] = f.SinusoidalImpulse(name='SinusMid', amplitude=1, phase_shift=0, wavelength=430.0e-09)
#test[20] = f.EnvelopeSinus(name='EnevelopedSinusMitte', wavelength=200e-09, phase_shift=0, amplitude=1, fwhm=600e-09, peak_timestep=200)
#test[50] = f.GaussianImpulse(name='a ndererImpuls', amplitude=1, peak_timestep=60, fwhm=100.e-09)
test[34] = f.ActivatedSinus(name='sin**2acivation', wavelength=900.0e-09, carrier_wavelength=6000.0e-9, phase_shift=0, amplitude=1)
test[96] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=6000)
test[0] = f.LeftSideMur()
test[100] = f.RightSideMur()

test.run_timesteps(6300)
test.get_observed_signals()


def theory_dielectric_slab(grid, d, n, timestep_start, timestep_end):
    k0 = grid.sources[0].omega / c0
    k = n * k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(-2j * k * d)
    theo_E = []
    theo_phase = []
    timestep_range = np.arange(timestep_start, timestep_end, 1)
    for ts in timestep_range:
        e_inc = np.exp(1j * grid.sources[0].omega * grid.dt * ts)
        e_tr = e_inc * (2/(n+1)) * (2*n / (n+1)) * (1/(1-q)) * np.exp(-1j*(k+k0)*d)
        amplitude = np.imag(e_tr) / 2
        theo_E.append(amplitude)
        theo_phase.append(np.angle(e_tr))
    return theo_E, theo_phase

data_theo_E, data_theo_phase = theory_dielectric_slab(grid=test, n=2, d=31*test.dz, timestep_start=6000, timestep_end=6300)


x = np.arange(6000, 6300, 1)
y = np.arange(0, 6300, 1)
z = np.arange(0, 6300, 1)

fig, axes = plt.subplots(1, 1)
axes.plot(x, data_theo_E, label='Class_theo')
axes.plot(z, test.local_observers[0].amplitude * np.cos(test.sources[0].omega * test.dt * y + test.local_observers[0].phase), label='reconstructed', linestyle='dashed')
axes.plot([test.local_observers[0].first_timestep, test.local_observers[0].second_timestep], test.local_observers[0].observedE, linestyle='None', marker='o')
axes.legend(loc='best')
plt.show()

print(data_theo_phase[0])