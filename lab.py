import fdtd_1d as f
import numpy as np

test = f.Grid(351, 4.3e-09) # creates 300 grid cells (á 6.3e-09m)

#test[80] = f.GaussianImpulse(name='GaussLinks', amplitude=1 ,peak_timestep=20, fwhm=80e-09)
#test[89] = f.GaussianImpulse(name='GaussLinks', amplitude=1 ,peak_timestep=20, fwhm=80e-09)
test[0] = f.LeftSideGridBoundary()
test[350] = f.RightSideGridBoundary()
#test[98] = f.RightSideGridBoundary()
#test[200:231] = f.NonDispersiveMedia(name='Medium4Eps', permeability=1, permittivity=4, conductivity=0)
test[200:271] = f.NonDispersiveMedia(name='Media2Epsleft', permeability=1, permittivity=2, conductivity=5e04)
#test[75:171] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=5e04)
#test[150:171] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=2, conductivity=0)
#test[20] = f.SinusoidalImpulse(name='SinusMid', amplitude=0.8, phase_shift=0, wavelength=630.0e-09)
#test[35] = f.EnvelopeSinus(name='EnevelopedSinusMitte', wavelength=100e-09, phase_shift=0, amplitude=1, fwhm=500e-09, peak_timestep=200)
test[30] = f.GaussianImpulse(name='andererImpuls', amplitude=1.3, peak_timestep=60, fwhm= 120.e-09)

test.animate_timesteps(1200)
