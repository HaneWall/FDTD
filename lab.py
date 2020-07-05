import fdtd_1d as f


'''
# Load setup/benchmark:
#setup = f.Harmonic_Slab_Setup(dx=4.e-09, length_grid_in_dx=40, length_media_in_dx=30, start_index_media=5, wavelength=240.e-09, epsilon=4, ampl=1, timesteps=2000)
#setup.run_benchmark()
#setup_2 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=25, wavelength_guided_for=800.e-09, wavelength=800e-09, ampl=1, timesteps=15000, number_of_layer_pairs=15, vary_layers=True, vary_inc_wavelength=False)
#setup_2.run_benchmark()
#setup_3 = f.Harmonic_Slab_Lorentz_Setup(dx=4/6 * 3.91e-07, length_grid_in_dx=50, length_media_in_dx=30, start_index_media=5, wavelength=1.75e-05, eps_inf=1.5, chi_1=[5.1], chi_2=[0], chi_3=[0], conductivity=0, w0=[1.2566e14], gamma=[8e12], ampl=1, timesteps=4000)
setup_3 = f.Harmonic_Slab_Lorentz_Setup(dx=3/6 * 3.91e-07, length_grid_in_dx=50, length_media_in_dx=30, start_index_media=5, wavelength=1.75e-05, eps_inf=1.2, chi_1=[2.1, 2.4], chi_2=[0, 0], chi_3=[0, 0], conductivity=0, w0=[1.2566e14, 1.2e13], gamma=[8e12, 9e14], ampl=1, timesteps=4000)
setup_3.run_benchmark()
'''

# Or build your own setup

# Step 1: init grid
test = f.Grid(100, 4e-09) # creates 201 grid cells (á 4.0e-09m)

# Step 2: init media
#test[30:80] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=3, conductivity=0)
#test[30:60] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[20:758] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[758:1496] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
test[20:50] = f.LorentzMedium(name='firsttest', permeability=1, eps_inf=1.5, chi_1=[2.1, 2.4], chi_2=[0, 0], chi_3=[0, 0], conductivity=0, w0=[1.2566e14, 1.2e15], gamma=[8e10, 9e10])
#test[180:190] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)


# Step 3: init sources
#test[15] = f.ActivatedSinus(name='sin**2activation', wavelength=200e-09, carrier_wavelength=1000.0e-09, phase_shift=0, amplitude=1, tfsf=True)
#test[10] = f.ActivatedSinus(name='sin**2activation', wavelength=1.4e-05, carrier_wavelength=20*1.75e-05, phase_shift=0, amplitude=1, tfsf=True)
#test[15] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=1000, tfsf=True)
#test[15] = f.EnvelopeSinus(name='test', wavelength=1.064e-06, fwhm=14.6e-06, amplitude=1, phase_shift=0, peak_timestep=3000, tfsf=True)
test[10] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=100, fwhm=60e-09, tfsf=True)
    
# Step 4: add observer
#test[280] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
#test[10] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
test[5] = f.FFTObserver(name='scndAttempt', first_timestep=0, second_timestep=10000)
test[80] = f.FFTObserver(name='firstattempt', first_timestep=0, second_timestep=10000)

# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[99] = f.RightSideMur()
    
# Step 6: run simulation
test.run_timesteps(timesteps=10000)
#test.animate_timesteps(60000)

# Step 7: misc
#print(test.dt)
#test.get_observed_signals()
#test.visualize_permittivity()
test.visualize_fft_observed()
#test.store_Ez_Data()