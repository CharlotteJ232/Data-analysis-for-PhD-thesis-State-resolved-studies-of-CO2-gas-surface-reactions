import numpy as np
from scipy import signal
from scipy.optimize import curve_fit


def energy_lockin(data, frequency, measure_freq, bandwidth, delta_t):
    """
    This function takes a voltage data array and calculates the energy at a
    given frequency and bandwidth for every time step delta_t.
    It subtracts the background energy by applying 2 extra lock-ins left
    and right of the cantilever energy. Since background subtraction is
    applied, negative energies can emerge.

    Returns the lock-in energy, sampled at 1/delta_t.

    Further documentation on lock-in amplifiers and used functions
    --------------------------------------------------------------
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/
    scipy.signal.butter.html

    https://en.wikipedia.org/wiki/Butterworth_filter

    https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/
    scipy.signal.filtfilt.html

    https://www.thinksrs.com/downloads/pdfs/applicationnotes/AboutLIAs.pdf

    Parameters
    ----------
    data : numpy 1D-array
        Voltage data from a .tdms file.

    frequency : float
        The fitted resonance frequency of the cantilever.

    measure_freq : float or int
        Frequency at which the data points in the .tdms file were measured.

    bandwidth : float or int
        Bandwidth for the lock-in at the resonance frequency of the cantilever.

    delta_t : float or int
        Time steps at which the cantilever energy is sampled. 1/100 Hz
        in the case of 100 Hz synchronization between temperature and .tdms
        files.

    Returns
    -------
    energy_min_b : numpy 1D-array
        1D-array with energy of the lock-in at the resonance frequency of the
        cantilever sampled at 1/delta_t, minus the mean background energy per
        Hz multiplied by the bandwidth.
    """
    data_points = len(data)  # the total number of data points

    sr = int(np.floor(delta_t * measure_freq))  # samples to read per 1/100 Hz

    # total energy points
    energy_points = int(np.floor(data_points / (delta_t * measure_freq)))

    # empty array for energy lock-in at the resonance frequency.
    energy_center = np.zeros(energy_points)

    # time array corresponding to the data points
    time_arr = np.linspace(0, data_points / measure_freq, data_points)

    # Calculate the energy at the res frequency ==============================
    # the cosine and sine at the frequency of the cantilever
    cosine = np.cos(2 * np.pi * frequency * time_arr)
    sine = np.sin(2 * np.pi * frequency * time_arr)

    # multiplying by cosine and sine gives a convolution with a delta in the
    # Fourier domain. The difference frequency will be zero, while the
    # sum frequency will lie outside the low pass filter range. So, effectively
    # half of the signal remains after the low pass filter.
    cos_data = data * cosine  # multiply data by cosine with cantilever freq
    sin_data = data * sine  # multiply data by sine with cantilever freq

    # cut-off frequency of a low pass butter filter. Bandwidth is twice the
    # cut-off frequency. Normalise by dividing by Nyquist freq.
    cutoff_freq = float(bandwidth / 2) / (measure_freq / 2)

    # find low pass digital Butterworth filter coefficients.
    b, a = signal.butter(2, cutoff_freq, 'low', analog=False, output='ba')

    # Apply low-pass filter to filter out other frequency components
    # lfilt is Python equivalent of Matlab filter function. Somehow the
    # Python function filtfilt does not produce the same results.
    # lfilter causes phase distortions, filtfilt does not.
    cosine_filter_data = signal.lfilter(b, a, cos_data)
    sine_filter_data = signal.lfilter(b, a, sin_data)

    # Calculate the energy at the resonance frequency of the cantilever
    # at every 1/delta_t step by taking the sum of the square of both.
    for i in range(energy_points):
        energy_center[i] = ((np.sum((
                            cosine_filter_data[i * sr:(i + 1) * sr] ** 2)
                            + sine_filter_data[i * sr:(i + 1) * sr] ** 2)
                            ) / measure_freq)

    # Subtract background: do lock-in 30 Hz to the left and right ============
    # empty array for the background energy
    energy_b = np.zeros((2, energy_points))
    index = 0  # row index in which the background energies will be stored

    # [Hz]. array with two shifts of lock-in frequency
    freq_shift = np.array([-30, 30])
    bandwidth_background = 20  # [Hz]. Take a larger box to average.

    # Determine cut-off frequency and low pass filter coefficients for
    # the background.
    cutoff_freq_b = (bandwidth_background / 2) / (measure_freq / 2)
    b, a = signal.butter(2, cutoff_freq_b, 'low', analog=False, output='ba')

    # for-loop for background energy lock-in 30 Hz to the left and the right.
    for shift in freq_shift:
        # [Hz]. Frequency for background lock-in.
        frequency_background = frequency + shift

        cosine_b = np.cos(2 * np.pi * frequency_background * time_arr)
        sine_b = np.sin(2 * np.pi * frequency_background * time_arr)

        # multiply data by frequency shifted lock-in
        cosine_data_b = data * cosine_b
        sine_data_b = data * sine_b

        # apply low pass filter to only inspect freqs near freq_background
        cosine_filter_b = signal.lfilter(b, a, cosine_data_b)
        sine_filter_b = signal.lfilter(b, a, sine_data_b)

        # Calculate the background energy at the given shifted freq and filter,
        # for every 1/delta_t step.
        for i in range(energy_points):
            energy_b[index, i] = ((np.sum(
                                    cosine_filter_b[i * sr:(i + 1) * sr] ** 2
                                    + sine_filter_b[i * sr:(i + 1) * sr] ** 2
                                    )) / measure_freq)

        index += 1  # add 1 to get the good row index for the next lock-in

    # Calculate average background energy per Hz (for 30 Hz shift lock-in).
    energy_b_left = np.mean(energy_b[0, :]) / bandwidth_background
    energy_b_right = np.mean(energy_b[1, :]) / bandwidth_background

    # Calculate the lock-in energy minus the average background energy
    energy_min_b = (energy_center
                    - np.mean([energy_b_left, energy_b_right]) * bandwidth)

    return energy_min_b


def exponential(xdata, p1, p2):
    # Exponential function for fit to the auto-correlation of the energy.
    return p1*np.exp(p2*xdata)


def get_qfactor(energy, res_freq, sampling_freq, filter_cutoff, q_guess=5e4):
    """
    Function that calculates the Q-factor and its estimated error using
    auto-correlation. For explanation see thesis Guide vd Stolpe, page 26/27.

    Parameters
    ----------
    energy : numpy 1D-array
        Lock-in energies, sampled at sampling_freq.

    res_freq : float or int
        Resonance frequency of the cantilever.

    sampling_freq : int or float
        Frequency at which the energy is sampled (100 Hz).

    filter_cutoff : int or float
        Bandwidth of the low pass filter used to calculate the lock-in energy
        at resonance frequency.

    q_guess : int or float
        Guess for the value of the Q-factor. Default is 50.000 .

    Returns
    -------
    q_factor : float or int
        Q-factor calculated with auto-correlation.

    q_factor_std : float or int
        Error estimation of calculated Q-factor.
    """
    max_time = len(energy) / sampling_freq  # [s]. Max time of energy array

    # [s]. Stop fitting to get estimated parameters. A longer time will give
    # a larger error. See thesis Guido vd Stolpe, page 26-28.
    stop_fit_time = 3.5

    # Due to the applied low pass filter (in the energy lock-in calculations)
    # we only need times of frequencies within the bandwidth.
    # Add factor 2 to be sure that the filter does not influence the result.
    filter_corr_time = 1./filter_cutoff

    # Subtract mean energy from energy array.
    relative_energy = energy - np.mean(energy)

    # Calculate auto correlation, using scipy. signal.fftconvolve because fft
    # method is way faster(!) than numpy.correlate etc. Use mode="full" to get
    # full correlation, see documentation of scipy.signal.fftconvolve.
    #
    # To use the fftconvolve function to calculate the cross correlation, we
    # must reverse the order of one of the arrays, because convolution differs
    # from correlation. So, reverse the second array.
    correlation = signal.fftconvolve(relative_energy, relative_energy[::-1],
                                     mode="full")

    # Take the one sided energy correlation function
    correlation_1 = correlation[round(len(correlation)/2)+1:]

    # Time lag corresponding to the auto-correlation function
    time_shift = np.linspace(0, len(energy)/sampling_freq, len(correlation_1))

    # List with start and stop index for fit of auto-correlation.
    fitrange = [int(round(sampling_freq*filter_corr_time)),
                int(round(sampling_freq*stop_fit_time))]

    # Magnify the data for fitting procedure with scipy.optimize_curvefit
    magnification_factor = 100./correlation_1[fitrange[0]]

    # Guess that the amplitude of the exponent is the maximum of the
    # one sided correlation multiplied by the magnification factor
    amplitude_guess = np.max(correlation_1*magnification_factor)

    p2_guess = -2. / q_guess  # see equation 3.6 and 3.7 in thesis Guido.
    param_guess = np.asarray([amplitude_guess, p2_guess])  # Guess for params.

    # Fit the auto-correlation of the energy. Calculate best fitting parameters
    # to exponential function, only using data from indices
    # [filter_corr_time*sampling_freq : stop_fit_time*sampling_freq]
    params, pcov = curve_fit(exponential,
                             time_shift[fitrange[0]:fitrange[1]],
                             magnification_factor
                             * correlation_1[fitrange[0]:fitrange[1]],
                             p0=param_guess)

    # Calculate characteristic time and corresponding Q-factor.
    tau = - 2 / params[1]
    q_factor = np.pi * res_freq * tau

    # "number of correlation times that fit into the cantilever energy file
    # as the correlation function is averaging over all these independent
    # samples." -- thesis Guido vd Stolpe, page 26.
    n_correlation_times = max_time / tau

    # Error estimation of Q-factor. See thesis Guido vd Stolpe, page 27.
    q_factor_std = 3.9 / np.sqrt(n_correlation_times) * q_factor

    return q_factor, q_factor_std
