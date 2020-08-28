import fdtd_1d as f
import numpy as np

# Load setup/benchmark:
#setup = f.Harmonic_Slab_Setup(dx=4.e-09, length_grid_in_dx=40, length_media_in_dx=30, start_index_media=5, wavelength=240.e-09, epsilon=4, ampl=1, timesteps=10000)
#setup.run_benchmark()
#setup_2 = f.TiO2_Si02_Dielectric_Mirror_Setup(N_lambda_media=25, wavelength_guided_for=800.e-09, wavelength=800e-09, ampl=1, timesteps=15000, number_of_layer_pairs=15, vary_layers=True, vary_inc_wavelength=False)
#setup_2.run_benchmark()
#setup_3 = f.Harmonic_Slab_Lorentz_Setup(name='new_numpy_step_P', dx=[4/6 * 3.91e-07], length_grid_in_dx=[50], length_media_in_dx=[30], start_index_media=5, wavelength=1.75e-05, eps_inf=1.5, chi_1=[5.1], chi_2=[0], chi_3=[0], conductivity=0, w0=[1.2566e14], gamma=[8e12], ampl=1, timesteps=[6000])
#setup_3.run_benchmark()

#setup_3 = f.Harmonic_Slab_Lorentz_Setup(name='new_data_files_per_terminal', dx=[2/6 * 3.91e-07, 4/6 * 3.91e-07, 6/6 * 3.91e-07], length_grid_in_dx=[100, 55, 40], length_media_in_dx=[90, 45, 30], start_index_media=5, wavelength=1.75e-05, eps_inf=1.05, chi_1=[2.1, 2.4], chi_2=[0, 0], chi_3=[0, 0], conductivity=0, w0=[1.2566e14, 1.2e13], gamma=[8e12, 9e14], ampl=1, timesteps=[15000, 7500, 5000], courant=1)
#setup_3.run_benchmark()
#setup_3.store_obs_data()

#setup5 = f.QPM_end_P(name='P_at_end_30000_1lambda', peak_timestep=24000, timesteps=60000, pulse_duration=20e-15, number_of_lambdas=1)
#setup5.run_benchmark()
#setup5.store_obs_data()

#setup6 = f.QPM_Length(name='1000obs_10fs_32000pt_05courant_6lambda', number_of_lambdas=6, timesteps=76000, peak_timestep=32000, pulse_duration=10e-15, number_of_distributed_observer=1000)
#setup6.run_benchmark()
#setup6.store_obs_data()
# Or build your own setup

# Step 1: init grid
test = f.Grid(3000, 4e-09, courant=0.5) # creates 201 grid cells (á 4.0e-09m)

# Step 2: init media
#test[20:40] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=3, conductivity=0)
#test[30:45] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
test[800:1000] = f.CentroRamanMedium(name='test', chi_1=[0.69617, 0.40794, 0.89748], w0=[2.7537e16, 1.6205e16, 1.9034e14], chi_3=[1.94e-22, 0, 0], alpha=[0.7, 0, 0], wr=[8.7722e13, 0, 0], gamma_K=[0, 0, 0], gamma_R=[3.1250e13, 0, 0], permeability=1, conductivity=0, eps_inf=1)

#test[800:1100] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[30:45] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[871:1722] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[1723:2574] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[2574:3425] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1.05, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])

#test[20:757] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[757:1494] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1, chi_1=[2.4272, 1.4617, 9.6536], chi_2=[-30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 7.9514e15, 9.7766e13], gamma=[0, 0, 0])
#test[1494:2231] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1, chi_1=[2.42, 9.65, 1.46], chi_2=[30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])
#test[2231:2968] = f.LorentzMedium(name='Varin', permeability=1, eps_inf=1, chi_1=[2.42, 9.65, 1.46], chi_2=[-30.e-12, 0, 0], chi_3=[0, 0, 0], conductivity=0, w0=[1.5494e16, 9.776e13, 7.9514e15], gamma=[0, 0, 0])


#test[20:50] = f.LorentzMedium(name='firsttest', permeability=1, eps_inf=1, chi_1=[2.1, 2.4], chi_2=[0, 0], chi_3=[0, 0], conductivity=0, w0=[1.2566e14, 1.2e15], gamma=[8e10, 9e10])
#test[180:190] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)
#test[15:25] = f.NonDispersiveMedia(name='randm', permeability=1, permittivity=1.4, conductivity=0)
#test[25:35] = f.NonDispersiveMedia(name='randm2', permeability=1, permittivity=1.4, conductivity=0)

# Step 3: init sources

#test[600] = f.SechEnveloped(name='first_try', wavelength=1.064e-06, pulse_duration=20e-15, Intensity=1, peak_timestep=60000, tfsf=False)
test[200] = f.GaussianImpulseWithFrequency(name='test', wavelength=1.5e-06, pulse_duration=10e-15, tfsf=False, Intensity=10e12, peak_timestep=16000)
#test[15] = f.ActivatedSinus(name='sin**2activation', wavelength=200e-09, carrier_wavelength=1000.0e-09, phase_shift=0, amplitude=1, tfsf=True)
#test[10] = f.ActivatedSinus(name='sin**2activation', wavelength=1.064e-06, carrier_wavelength=20*1.064e-06, phase_shift=0, amplitude=1, tfsf=True)
#test[5] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=80.0e-09, tfsf=True)

#test[600] = f.EnvelopeSinus(name='test', wavelength=1.064e-06, fwhm=14.6e-06, amplitude=2*2.74492e7, phase_shift=0, peak_timestep=16000, tfsf=False)

#test[3] = f.EnvelopeSinus(name='test', wavelength=1.064e-06, fwhm=14.6e-06, amplitude=2*2.74492e7, phase_shift=0, peak_timestep=9000, tfsf=True)
#test[20] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=100, fwhm=30e-09, tfsf=False)
#test[20] = f.GaussianImpulseWithFrequency(name='varin', Intensity=10e08, wavelength=1.064e-06, pulse_duration=20e-15, peak_timestep=20000, tfsf=True)
#test[5] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=5000, fwhm=30e-09/5, tfsf=True)

# Step 4: add observer

#test[280] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
#test[10] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
#test[5] = f.E_FFTObserver(name='scndAttempt', first_timestep=0, second_timestep=10000)
#test[50] = f.E_FFTObserver(name='just_testing_2', first_timestep=0, second_timestep=22000)

#test[2972] = f.E_FFTObserver(name='E_two_lambda_laterpeak_16000_737', first_timestep=0, second_timestep=99999)
#test[2966] = f.P_FFTObserver(name='P_two_lambda_laterpeak_16000_737', first_timestep=0, second_timestep=99999)

# Step 5: add boundaries
#test[0] = f.LeftSideMur()
test[0] = f.LeftSideGridBoundary()
test[2999] = f.RightSideGridBoundary()
#test[2999] = f.RightSideMur()

#test[3449] = f.RightSideMur()
#test[1499] = f.RightSideMur()
    
# Step 6: run simulation
test.run_timesteps(20000)
#test.animate_timesteps(60000)

# Step 7: misc
#test.get_observed_signals()
#test.visualize_permittivity()
#test.visualize_fft_observed()
#test.store_obs_data()
