import numpy as np

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

