import fdtd_1d as f
import matplotlib.pyplot as plt
import numpy as np

# Load setup/benchmark:
#setup = f.Harmonic_Slab_Setup(dx=4.e-09, length_grid_in_dx=40, length_media_in_dx=30, start_index_media=5, wavelength=160.e-09, epsilon=4, ampl=1, timesteps=2000)
#setup.run_benchmark()
#setup_2 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=40, wavelength_guided_for=800.e-09, wavelength=800e-09, ampl=1, timesteps=5500, number_of_layer_pairs=5, vary_layers=False, vary_inc_wavelength=False)
#setup_2.run_benchmark()


incident_wavelengths = [400.e-09 + (i / 201) * 800.e-09 for i in np.arange(0, 201, 1)]
setup_2 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=10, wavelength_guided_for=800.e-09, wavelength=400e-09, ampl=1, timesteps=3000, number_of_layer_pairs=5, vary_layers=False, vary_inc_wavelength=True)
setup_2.run_benchmark()

dbr_2 = np.array(setup_2.refl_ampl)**2
setup_5 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=10, wavelength_guided_for=800.e-09, wavelength=400e-09, ampl=1, timesteps=4000, number_of_layer_pairs=5, vary_layers=False, vary_inc_wavelength=True)
setup_5.run_benchmark()
dbr_5 = np.array(setup_5.refl_ampl)**2
setup_10 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=10, wavelength_guided_for=800.e-09, wavelength=400e-09, ampl=1, timesteps=4000, number_of_layer_pairs=10, vary_layers=False, vary_inc_wavelength=True)
setup_10.run_benchmark()
dbr_10 = np.array(setup_10.refl_ampl)**2

fig, axes = plt.subplots(1, 1)
axes.plot(np.array(incident_wavelengths)*10**9, dbr_2, linestyle='dashed', color='red', marker='o', alpha=0.5, label='DBR-2')
axes.plot(np.array(incident_wavelengths)*10**9, dbr_5, linestyle='dashed', color='blue', marker='o', alpha=0.5, label='DBR-5')
axes.plot(np.array(incident_wavelengths)*10**9, dbr_10, linestyle='dashed', color='green', marker='o', alpha=0.5, label='DBR-10')
axes.legend(loc='best')
axes.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
axes.set_xlabel(r'$\lambda$ in nm', fontsize=14)
axes.set_ylabel(r'Reflektionsgrad $\mathcal{R}$', fontsize=14)
plt.show()


'''
refl_amplitudes = []
theorie = []
layer_numbers = range(1, 21)

for layer in layer_numbers:
    setup = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=25, wavelength_guided_for=800.e-09, wavelength=800e-09, ampl=1, timesteps=7000, number_of_layer_pairs=layer)
    setup.run_benchmark()
    refl_amplitudes.append(setup.refl_ampl)
    theorie.append(((setup.si_n**(2*layer)-setup.ti_n**(2*layer))/(setup.si_n**(2*layer)+setup.ti_n**(2*layer)))**2)

fig, axes = plt.subplots(1, 1)
axes.plot(layer_numbers, np.array(refl_amplitudes)**2, linestyle='dashed', color='blue', marker='o', alpha=0.5, label='FDTD')
axes.plot(layer_numbers, np.array(theorie), linestyle='dashed', color='red', marker='o', alpha=0.5, label='model')
axes.legend(loc='best')
axes.grid(True, linestyle=(0, (1, 5)), color='black', linewidth=1)
axes.set_xlabel('Anzahl der Paare', fontsize=14)
axes.set_ylabel(r'Reflektionsgrad $\mathcal{R}$', fontsize=14)
plt.show()


# Or build your own setup

# Step 1: init grid
test = f.Grid(101, 4.0e-09) # creates 201 grid cells (á 4.0e-09m)

# Step 2: init media
test[20:51] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)
test[51:90] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)

# Step 3: init sources
#test[15] = f.ActivatedSinus(name='sin**2activation', wavelength=80.0e-09, carrier_wavelength=10000.0e-9, phase_shift=0, amplitude=1, tfsf=True)
#test[49] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=80.e-09, tfsf=True)
#test[10] = f.EnvelopeSinus(name='test', wavelength=60.0e-09, fwhm=200.e-09, amplitude=1, phase_shift=0, peak_timestep=300,tfsf=True)
test[15] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=120, fwhm=130.e-09, tfsf=True)
    
# Step 4: add observer
#test[90] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
    
# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[100] = f.RightSideMur()
    
# Step 6: run simulation
#test.run_timesteps(timesteps=2200, vis=False)
test.animate_timesteps(2000)

# Step 7: misc
'''