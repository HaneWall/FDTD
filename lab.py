import fdtd_1d as f
import numpy as np
#import matplotlib
#matplotlib.use('qt4agg')
import matplotlib.pyplot as plt

from fdtd_1d.constants import c0


def theory_dielectric_slab_new(grid):
    d = ((grid.materials[0].position[-1] - grid.materials[0].position[0] + 1)) * grid.dx
    n = np.sqrt(grid.eps[grid.materials[0].position[1]])
    k0 = grid.sources[0].omega / c0
    k = n*k0
    q = ((n - 1) ** 2) / ((n + 1) ** 2) * np.exp(2j * k * d)
    e_inc = grid.sources[0].ampl * np.exp(-1j*np.pi/2)
    e_tr = e_inc * (2/(n+1)) * (2*n / (n+1)) * (1/(1-q)) * np.exp(1j*(k-k0)*d)
    theo_amplitude = np.abs(e_tr)
    theo_phasenunterschied = np.angle(e_tr)
    return theo_amplitude, theo_phasenunterschied

# BUILD SETUP
indices = [52 + i for i in np.arange(0, 20, 1)]
theo_phasenunterschied = []
theo_amplitude = []
exp_phase = []
exp_amplitude = []

for ind in indices:
    # Step 1: init grid
    test = 'test'+str(ind)
    test = f.Grid(201, 4.0e-09) # creates 201 grid cells (á 4.3e-09m)

    # Step 2: init media
    test[50:ind] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)

    # Step 3: init sources
    #test[49] = f.ActivatedSinus(name='sin**2activation', wavelength=113.6e-09, carrier_wavelength=10000.0e-9, phase_shift=0, amplitude=1, tfsf=True)
    test[49] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=100.e-09, tfsf=True)
    #test[20] = f.EnvelopeSinus(name='test', wavelength=200.0e-09, fwhm=600.e-09, amplitude=1, phase_shift=0, peak_timestep=300)
    #test[20] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=120, fwhm=130.e-09)
    
    # Step 4: add observer
    test[190] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
    
    # Step 5: add boundaries
    test[0] = f.LeftSideMur()
    test[200] = f.RightSideMur()
    
    # Step 6: run simulation
    test.run_timesteps(2200)
    
    # Step 7: misc
    exp_amplitude.append(test.local_observers[0].amplitude)
    exp_phase.append(test.local_observers[0].phase)
    theo_amplitude.append(theory_dielectric_slab_new(test)[0])
    theo_phasenunterschied.append(theory_dielectric_slab_new(test)[1])


fig, axes = plt.subplots(2, 2)
fig.suptitle(r'$N_{\lambda_{vac}}=10$', fontsize=20)
axes[0][0].plot(np.array(indices)-51, np.array(theo_amplitude), label='theorie', color='blue')
axes[0][0].grid(True, linestyle='-', alpha=0.4)
axes[0][0].plot(np.array(indices)-51, np.array(exp_amplitude), label='gemessen', linestyle='dashed', color='red')
axes[0][0].legend(loc='best')
axes[0][0].set_xlabel('Breite des Mediums in '+r'$\Delta_x$')
axes[0][0].set_ylabel('Transmittierte Amplitude '+ r'$Ez_{tr}$')
axes[0][1].plot(np.array(indices)-51, np.array(theo_amplitude)/np.array(exp_amplitude), color='black')
axes[0][1].set_ylabel(r'$\frac{E_{tr,theo}}{E_{tr,exp}}$')
axes[0][1].set_xlabel('Breite des Mediums in '+r'$\Delta_x$')
axes[0][1].grid(True, linestyle='-', alpha=0.4)
axes[1][0].set_ylabel('Phasenverlauf')
axes[1][0].plot(np.array(indices)-51, theo_phasenunterschied, label='Theoretischer Phasenuntserschied', color='blue')
axes[1][0].set_xlabel('Breite des Mediums in '+r'$\Delta_x$')
axes[1][0].grid(True, linestyle='-', alpha=0.4)
axes[1][0].plot(np.array(indices)-51, -np.array(exp_phase), label='gemessene Phase', color='red', linestyle='dashed')
axes[1][0].legend()
axes[1][1].set_xlabel('Breite des Mediums in '+r'$\Delta_x$')
axes[1][1].set_ylabel(r'$d(\phi_{exp},\phi_{theo})$')
axes[1][1].plot(np.array(indices)-50, np.abs(-np.array(exp_phase) - np.array(theo_phasenunterschied)), color='black')
axes[1][1].grid(True, linestyle='-', alpha=0.4)
plt.show()

'''
# Step 1: init grid
test = f.Grid(nx=201, dx=4.0e-09) # creates 201 grid cells (á 4.3e-09m)

# Step 2: init media
test[100:101] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)

# Step 3: init sources
test[60] = f.ActivatedSinus(name='sin**2acivation', wavelength=300.0e-09, carrier_wavelength=10000.0e-9, phase_shift=0, amplitude=1, tfsf=True)
#test[20] = f.EnvelopeSinus(name='test', wavelength=80.0e-09, fwhm=150.e-09, amplitude=1, phase_shift=0, peak_timestep=300, tfsf=True)
#test[30] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=300, fwhm=200.e-09, tfsf=True)
#test[20] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=70.e-09, tfsf=True)
# Step 4: add observer

#test[190] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=800)

# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[200] = f.RightSideMur()

# Step 6: run simulation
#test.run_timesteps(900)
# Step 7: Misc
#test.get_observed_signals()

print(len(test.eps))
'''