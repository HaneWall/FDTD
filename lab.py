import fdtd_1d as f
import numpy as np

# Load setup/benchmark:
setup = f.Harmonic_Slab_Setup(dx=4.e-09, length_grid_in_dx=50, length_media_in_dx=30, start_index_media=5, wavelength=90.e-09, epsilon=4, ampl=1, timesteps=1800)
setup.run_benchmark()

'''
# Or build your own setup

# Step 1: init grid
test = f.Grid(101, 4.0e-09) # creates 201 grid cells (á 4.0e-09m)

# Step 2: init media
test[20:81] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)

# Step 3: init sources
test[15] = f.ActivatedSinus(name='sin**2activation', wavelength=80.0e-09, carrier_wavelength=10000.0e-9, phase_shift=0, amplitude=1, tfsf=False)
#test[49] = f.SinusoidalImpulse(name='test', amplitude=1, phase_shift=0, wavelength=80.e-09, tfsf=True)
#test[20] = f.EnvelopeSinus(name='test', wavelength=200.0e-09, fwhm=600.e-09, amplitude=1, phase_shift=0, peak_timestep=300)
#test[20] = f.GaussianImpulse(name='test', amplitude=1, peak_timestep=120, fwhm=130.e-09)
    
# Step 4: add observer
test[90] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=2000)
    
# Step 5: add boundaries
test[0] = f.LeftSideMur()
test[100] = f.RightSideMur()
    
# Step 6: run simulation
#test.run_timesteps(timesteps=2200, vis=False)
test.animate_timesteps(2000)

# Step 7: misc
'''