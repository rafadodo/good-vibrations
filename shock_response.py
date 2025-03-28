# -*- coding: utf-8 -*-
"""
This module contains functions for the computation of shock responses and response spectra.

Functions:
- get_shock_response_spectrum_from_signal(...): Obtains the shock response spectrum for a given
    input acceleration history, natural frequency range, and damping ratio.
- get_shock_response(...): Computes the acceleration response for a given single-degree-of-freedom
system and acceleration input.
"""
import numpy as np
import scipy as sp


def get_shock_response_spectrum(input_accel_history,
                                sampling_freq,
                                damping_ratio,
                                freq_range,
                                freq_resolution):
    """Obtains the shock response spectrum that corresponds to a given input acceleration time
    history, by computing the maximum response it produces in a range of single-degree-of-freedom
    systems, defined by a shared damping ratio and a range of natural frequencies. The method used
    is the filter scheme proposed by Smallwood as described in get_shock_response().

    Parameters:
        input_accel_history (numpy.ndarray): Input acceleration time history [g].
        sampling_freq (float): Sampling frequency [Hz] of the 
        damping_ratio (float): Node numbers to search for in the file.
        freq_range (list): Upper and lower boundaries of the output spectrum frequency range [Hz].
        freq_resolution (float): Frequency resolution [Hz] of the output spectrum.

    Returns:
        frqs (numpy.ndarray): Frequency array for the response spectrum [Hz].
        response_spectrum (numpy.ndarray): Shock response spectrum [g].
    """
    freqs = np.arange(freq_range[0], freq_range[1], freq_resolution)
    response_spectrum = np.zeros(len(freqs))
    for index in range(len(freqs)):
        response_spectrum[index] = max(abs(get_shock_response(input_accel_history,
                                                            sampling_freq,
                                                            freqs[index],
                                                            damping_ratio)))

    return freqs, response_spectrum


def get_shock_response(input_accel_history, sampling_freq, natural_freq, damping_ratio):
    """Computes the response acceleration time history for a single-degree-of-freedom system,
    defined by its natural frequency and damping ratio, when excited by a given input acceleration
    time history. The method used is the ramp-invariant recursive filtering proposed by Smallwood
    (see Reference).

    Parameters:
        input_accel_history (numpy.ndarray): Input acceleration time history [g].
        sampling_freq (float): Sampling frequency [Hz] of the input time history.
        natural_freq (float): Natural frequency [Hz] of the 1DOF system to evaluate.
        damping_ratio (float): Damping ratio of the 1DOF system to evaluate.
        output spectrum must be delivered.

    Returns:
        numpy.ndarray: Response acceleration time history [g].

    Reference:
        Smallwood, David O. “Improved Recursive Formula for Calculating Shock Response Spectra”
        Shock and Vibration Bulletin 51, no. 2, 1980.
    """
    delta_t = 1 / sampling_freq # sampling interval [s]
    omega_n = 2 * np.pi * natural_freq # natural angular frequency (rad/s)
    zeta = damping_ratio
    omega_d = omega_n*np.sqrt(1-zeta**2) # damped natural angular frequency (rad/s)

    E = np.exp(-zeta * omega_n * delta_t)
    C = E * np.cos(omega_d * delta_t)
    S = E * np.sin(omega_d * delta_t)
    S_d = S / omega_d / delta_t

    b = np.zeros(3) # filter numerator coefficients array
    a = np.zeros(3) # filter denominator coefficients array
    b[0]= 1 - S_d
    b[1] = 2* (S_d - C)
    b[2] = E**2 - S_d
    a[0] = 1
    a[1] = -2 * C
    a[2] = E**2

    response_accel_history = sp.signal.lfilter(b, a, input_accel_history)

    return response_accel_history
