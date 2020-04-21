import fdtd_1d as f
import numpy as np
import matplotlib.pyplot as plt

from fdtd_1d.constants import c0, mu0
from fdtd_1d.utilities import get_amplitude_and_phase

'''
def theory_dielectric_slab_new(grid):
    d = ((grid.materials[0].position[-1] - grid.materials[0].position[0]) + 1) * grid.dx
    n = np.sqrt(grid.eps[grid.materials[0].position[1]])
    k0 = grid.sources[0].omega / c0
    k = n*k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(2j * k * d)
    e_inc = grid.sources[0].ampl/2
    #e_inc = 0.506233
    e_tr = e_inc * (2/(n+1)) * (2*n / (n+1)) * (1/(1-q)) * np.exp(1j*(k-k0)*d)
    theo_amplitude = np.abs(e_tr)
    theo_phasenunterschied = np.angle(e_tr)
    return theo_amplitude, theo_phasenunterschied

# BUILD SETUP
indices = [81 + i for i in np.arange(0, 20, 1)]
theo_phasenunterschied = []
theo_amplitude = []
exp_phase = []
exp_amplitude = []

for ind in indices:
    # Step 1: init grid
    test = 'test'+str(ind)
    test = f.Grid(201, 4.0e-09) # creates 201 grid cells (á 4.3e-09m)

    # Step 2: init media
    test[50:ind] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)

    # Step 3: init sources
    test[40] = f.ActivatedSinus(name='sin**2acivation', wavelength=60.0e-09, carrier_wavelength=10000.0e-9, phase_shift=0, amplitude=1)
    #test[20] = f.EnvelopeSinus(name='test', wavelength=200.0e-09, fwhm=600.e-09, amplitude=1, phase_shift=0, peak_timestep=300)
    #test[20] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=120, fwhm=130.e-09)
    # Step 4: add observer
    test[190] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=1400)

    # Step 5: add boundaries
    test[0] = f.LeftSideMur()
    test[200] = f.RightSideMur()
    # Step 6: run simulation
    test.run_timesteps(1600)
    # Step 7: misc
    exp_amplitude.append(test.local_observers[0].amplitude)
    exp_phase.append(-test.local_observers[0].phase)
    theo_amplitude.append(theory_dielectric_slab_new(test)[0])
    theo_phasenunterschied.append(theory_dielectric_slab_new(test)[1])


fig, axes = plt.subplots(3, 1)
axes[0].plot(np.array(indices)-50, theo_amplitude, label='Theo')
axes[0].plot(np.array(indices)-50, exp_amplitude)
axes[0].legend()
axes[1].plot(np.array(indices)-50, theo_phasenunterschied, label='Theo')
axes[1].plot(np.array(indices)-50, np.array(exp_phase))
axes[1].legend()
axes[2].plot(np.array(indices)-50, np.array(exp_phase) - np.array(theo_phasenunterschied))
plt.show()

'''
# Step 1: init grid
test = f.Grid(nx=201, dx=4.0e-09) # creates 201 grid cells (á 4.3e-09m)

# Step 2: init media
#test[50:82] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)

# Step 3: init sources
test[20] = f.ActivatedSinus(name='sin**2acivation', wavelength=60.0e-09, carrier_wavelength=1000.0e-9, phase_shift=0, amplitude=1)
#test[20] = f.EnvelopeSinus(name='test', wavelength=200.0e-09, fwhm=600.e-09, amplitude=1, phase_shift=0, peak_timestep=300)
#test[30] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=120, fwhm=100.e-09)
# Step 4: add observer
test[90] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=800)

# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[200] = f.RightSideMur()
# Step 6: run simulation
test.run_timesteps(1000)
# Step 7: Misc
test.get_observed_signals()

