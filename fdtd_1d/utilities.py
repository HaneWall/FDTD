import numpy as np
from fdtd_1d.constants import c0

# get phase and amplitude based on Acos(omega * t + phase) with A > 0
def get_amplitude_and_phase(grid, first_timestep, second_timestep, data):
    signed_phase = np.arctan2(data[1] * np.cos(grid.sources[0].omega * grid.dt * first_timestep) - data[0] * np.cos(grid.sources[0].omega * grid.dt * second_timestep),
                      data[1] * np.sin(grid.sources[0].omega * grid.dt * first_timestep) - data[0] * np.sin(grid.sources[0].omega * grid.dt * second_timestep))

    signed_amplitude = data[0] / (
                np.cos(grid.sources[0].omega * grid.dt * first_timestep) * np.cos(signed_phase) -
                np.sin(grid.sources[0].omega * grid.dt * first_timestep) * np.sin(signed_phase))

    if signed_amplitude >= 0:
        amplitude = signed_amplitude
        phase = signed_phase
        return amplitude, phase

    elif signed_amplitude < 0 and signed_phase >= 0:
        amplitude = -signed_amplitude
        phase = signed_phase - np.pi
        return amplitude, phase

    elif signed_amplitude < 0 and signed_phase < 0:
        amplitude = -signed_amplitude
        phase = signed_phase + np.pi
        return amplitude, phase

def numerical_group_velocity(central_wavelength, dx, courant, group_velocity, n_real):
    omega = 2*np.pi*c0/central_wavelength
    numerical_g_v = group_velocity * (np.cos((n_real * omega * dx)/(2 * c0))) / (np.sqrt(1 - (group_velocity*courant/c0)**2 * np.sin(n_real * omega * dx /(2*c0))**2))
    return numerical_g_v

'''group_velocity = 205015332.2021728
central_wavelength = 1.5e-6
dx = 25e-9
courant = 0.5
n_real = 1.4446181264931932

print(numerical_group_velocity(central_wavelength=central_wavelength, dx=dx, courant=courant, group_velocity=group_velocity, n_real=n_real)/group_velocity)'''