import fdtd_1d as f

#   Spielwiese

test = f.Grid(201, 4.3e-09) # creates 201 grid cells (á 4.3e-09m)

test[0] = f.LeftSideGridBoundary()
test[200] = f.RightSideGridBoundary()   # note that 201 grid cells means cell 200 is the last
#test[80:121] = f.NonDispersiveMedia(name='Medium4Eps', permeability=1, permittivity=4, conductivity=0)
#test[100:171] = f.NonDispersiveMedia(name='Media2Epsleft', permeability=1, permittivity=2, conductivity=5e04)
#test[75:171] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=5e04)
test[90:111] = f.NonDispersiveMedia(name='Media4Epsleft', permeability=1, permittivity=4, conductivity=0)
#test[20] = f.SinusoidalImpulse(name='SinusMid', amplitude=0.8, phase_shift=0, wavelength=630.0e-09)
#test[35] = f.EnvelopeSinus(name='EnevelopedSinusMitte', wavelength=100e-09, phase_shift=0, amplitude=1, fwhm=500e-09, peak_timestep=200)
#test[30] = f.GaussianImpulse(name='andererImpuls', amplitude=1, peak_timestep=60, fwhm=120.e-09)
test[20] = f.ActivatedSinus(name='sin**2acivation', wavelength=230.0e-09, carrier_wavelength= 3000.0e-9, phase_shift=0, amplitude=1)
test[180] = f.QuasiHarmonicObserver(name='firstobserver', first_timestep=800)
test.animate_timesteps(1200)
print(test.local_observers[0].second_timestep)

