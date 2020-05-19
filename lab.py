import fdtd_1d as f

# Load setup/benchmark:
#setup = f.Harmonic_Slab_Setup(dx=4.e-09, length_grid_in_dx=40, length_media_in_dx=30, start_index_media=5, wavelength=320.e-09, epsilon=4, ampl=1, timesteps=2000)
#setup.run_benchmark()
#setup_2 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=20, wavelength_guided_for=800.e-09, wavelength=800e-09, ampl=1, timesteps=10000, number_of_layer_pairs=10, vary_layers=True, vary_inc_wavelength=False)
#setup_2.run_benchmark()


# Or build your own setup

# Step 1: init grid
test = f.Grid(201, 5e-07) # creates 201 grid cells (á 4.0e-09m)

# Step 2: init media
#test[50:180] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)
test[100:171] = f.LorentzMedium(name='firsttest', permeability=1, eps_inf=10, gamma=5e12, w0=100e12, chi_1=2.1, conductivity=0)
#test[51:90] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)

# Step 3: init sources
#test[15] = f.ActivatedSinus(name='sin**2activation', wavelength=200e-09, carrier_wavelength=1000.0e-09, phase_shift=0, amplitude=1, tfsf=True)
test[15] = f.ActivatedSinus(name='sin**2activation', wavelength=2.05e-05, carrier_wavelength=50e-05, phase_shift=0, amplitude=1, tfsf=True)
#test[15] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=1000, tfsf=True)
#test[15] = f.EnvelopeSinus(name='test', wavelength=1.9e-05, fwhm=3e-05, amplitude=1, phase_shift=0, peak_timestep=200,tfsf=True)
#test[15] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=300, fwhm=600, tfsf=True)
    
# Step 4: add observer
#test[90] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
#test[10] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[200] = f.RightSideMur()
    
# Step 6: run simulation
#test.run_timesteps(timesteps=2200, vis=False)
test.animate_timesteps(2000)

# Step 7: misc
#test.get_observed_signals()

