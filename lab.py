import fdtd_1d as f
import numpy as np
import matplotlib.pyplot as plt

from fdtd_1d.constants import c0
from fdtd_1d.utilities import get_amplitude_and_phase


def theory_dielectric_slab(grid, d, n, timestep_start, timestep_end):
    k0 = grid.sources[0].omega / c0
    k = n * k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(-2j * k * d)
    theo_E = []
    timestep_range = np.arange(timestep_start, timestep_end, 1)
    for ts in timestep_range:
        e_inc = np.exp(1j * grid.sources[0].omega * grid.dt * ts)
        e_tr = e_inc * (2/(n+1)) * (2*n / (n+1)) * (1/(1-q)) * np.exp(-1j*(k+k0)*d)
        amplitude = np.imag(e_tr) / 2
        theo_E.append(amplitude)
        theo_amplitude, theo_phase = get_amplitude_and_phase(grid=grid, first_timestep=timestep_start, second_timestep=(timestep_end - 1), data=[theo_E[0], theo_E[-1]])
    return theo_E, theo_amplitude, theo_phase


# BUILD SETUP

# Step 1: init grid
test = f.Grid(301, 4.3e-09) # creates 201 grid cells (á 4.3e-09m)

# Step 2: init media
test[200:251] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)

# Step 3: init sources
#test[20] = f.ActivatedSinus(name='sin**2acivation', wavelength=900.0e-09, carrier_wavelength=6000.0e-9, phase_shift=0, amplitude=1)
test[20] = f.EnvelopeSinus(name='test', wavelength=200.0e-09, fwhm=600.e-09, amplitude=1, phase_shift=0, peak_timestep=300)
# Step 4: add observer
#test[96] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=6000)

# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[300] = f.RightSideMur()

# Step 6: run simulation
test.animate_timesteps(2000)

# Step 7: misc
#test.get_observed_signals()
'''

wavelengths = [100.e-09 + i * 50.e-09 for i in range(0, 31)]
measured_phase = []
measured_amplitude = []
theo_phase = []
theo_amplitude = []
N_lambda = []

for wavelength in wavelengths:
    obj = 'wave' + str(wavelength)
    obj = f.Grid(101, 4.3e-09)
    obj[50:81] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)
    obj[34] = f.ActivatedSinus(name='sin**2acivation', wavelength=wavelength, carrier_wavelength=6000.0e-9, phase_shift=0, amplitude=1)
    obj[96] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=6000)
    obj[0] = f.LeftSideMur()
    obj[100] = f.RightSideMur()
    obj.run_timesteps(8000)
    measured_phase.append(obj.local_observers[0].phase)
    measured_amplitude.append(obj.local_observers[0].amplitude)
    theo_amplitude.append(theory_dielectric_slab(grid=obj, n=2, d=31*obj.dz, timestep_start=6000, timestep_end=6300)[1])
    theo_phase.append(theory_dielectric_slab(grid=obj, n=2, d=31*obj.dz, timestep_start=6000, timestep_end=6300)[2])
    N_lambda.append(wavelength/obj.dz)

fig, axes = plt.subplots(1, 2)

axes[0].plot(N_lambda, np.array(measured_amplitude)/np.array(theo_amplitude), marker='o', linestyle = 'dashed', label='measured/theoretical amplitude')
axes[0].plot(N_lambda, np.ones(len(wavelengths)), label='exact')
axes[0].legend(loc='best')
axes[0].set_xlabel(r'$N_{\lambda}$')
axes[1].plot(N_lambda, np.array(measured_phase)/np.array(theo_phase), marker='o', linestyle='dashed', label='measured/theoretical phase')
axes[1].plot(N_lambda, np.ones(len(wavelengths)), label='exact')
axes[1].legend(loc='best')
axes[1].set_xlabel(r'$N_{\lambda}$')
plt.show()


data_theo_E, data_theo_amplitude, data_theo_phase = theory_dielectric_slab(grid=test, n=2, d=31*test.dz, timestep_start=6000, timestep_end=6300)


x = np.arange(6000, 6300, 1)
y = np.arange(0, 6300, 1)
z = np.arange(0, 6300, 1)

fig, axes = plt.subplots(1, 1)
axes.plot(x, data_theo_E, label='Class_theo')
axes.plot(z, test.local_observers[0].amplitude * np.cos(test.sources[0].omega * test.dt * y + test.local_observers[0].phase), label='reconstructed', linestyle='dashed')
axes.plot([test.local_observers[0].first_timestep, test.local_observers[0].second_timestep], test.local_observers[0].observedE, linestyle='None', marker='o')
axes.legend(loc='best')
plt.show()

print(data_theo_amplitude, data_theo_phase)
'''