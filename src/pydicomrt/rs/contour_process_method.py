# import numpy as np
# import math

# from scipy import fft

# def contour_process(x_points, y_points, hierarchy, *args, **kwargs):
#     # print(x_points, y_points, hierarchy, args, kwargs)
#     if len(x_points) == 0 or len(y_points) == 0: return x_points, y_points  # check if contour is empty
    
#     if ("external_noise_size" in kwargs) and hierarchy == -1:
#         external_noise_size = kwargs["external_noise_size"]
#         x_points, y_points = ctr_external_denoise(x_points, y_points, noise_size=external_noise_size)
#         if len(x_points) == 0 or len(y_points) == 0: return x_points, y_points

#     if ("internal_noise_size" in kwargs) and hierarchy != -1:
#         internal_noise_size = kwargs["internal_noise_size"]
#         x_points, y_points = ctr_internal_denoise(x_points, y_points, noise_size=internal_noise_size)
#         if len(x_points) == 0 or len(y_points) == 0: return x_points, y_points
    
#     if "low_pass_ratio" in kwargs:
#         low_pass_ratio = kwargs["low_pass_ratio"]
#         x_points, y_points = ctr_fft_smooth(x_points, y_points, low_pass_ratio)

#     return x_points, y_points


# def ctr_internal_denoise(x_points, y_points, noise_size):
#     if len(x_points) < noise_size or len(y_points) < noise_size:
#         x_points = []
#         y_points = []
#     return x_points, y_points


# def ctr_external_denoise(x_points, y_points, noise_size):
#     if len(x_points) < noise_size or len(y_points) < noise_size:
#         x_points = []
#         y_points = []
#     return x_points, y_points


# def ctr_fft_smooth(x_points, y_points, low_pass_ratio=10):
#     contour_len = len(x_points)
#     if contour_len < 3:             # check if contour is too short
#         return x_points, y_points
    
#     x_points, y_points = resample_2d_contour(x_points, y_points)
#     x_points, y_points = fft_low_pass(x_points, y_points, filter_ratio=low_pass_ratio)
#     x_points = x_points[:contour_len]
#     y_points = y_points[:contour_len]
#     return x_points, y_points


# def resample_2d_contour(x, y, min_pt=1024):
#     repeat = math.ceil(2 + (min_pt / len(x)))
#     new_x = np.tile(x, repeat)
#     new_y = np.tile(y, repeat)
#     return new_x, new_y


# def fft_low_pass(x, y, filter_ratio=8):
#     """Apply low pass filter to x and y arrays by zeroing out high freq components.
    
#     ** actually, it's a band pass filter **

#     Parameters
#     ----------
#     x, y : array_like
#         Input arrays.
#     length : int > len(x)
#         The pass frequency length in frequency domain

#     Returns
#     -------
#     np.real(newZ), np.imag(newZ) : array_like
#         Filtered arrays of x and y
#     """
#     if filter_ratio >= 45:
#         return x, y

#     z = x + 1j * y    # Let x, y info together in signal space
#     C = fft.fft(z)  # R -> C

#     demod = np.abs(C)
#     indices = np.where(demod > 10)[0]

#     wave_length = len(x)
#     pass_length = int(wave_length * filter_ratio * 0.01)
#     mask = np.zeros_like(C, dtype=float)
#     length = pass_length + 1

#     # fix the bug if filtered signal's main freq only 0
#     try:
#         if length < indices[1]:
#             length = indices[1] + 1
#     except IndexError:
#         print("no other freq. > 10")
#         length = int(wave_length * 10 * 0.01)
    
#     mask[0] = 1
#     mask[1:length] = 1
#     mask[-length:] = 1

#     C *= mask   # apply filter mask
    
#     newZ = fft.ifft(C)  # C -> R
#     return np.real(newZ), np.imag(newZ)

import numpy as np
import math

def contour_process(x_points, y_points, hierarchy, *args, **kwargs):
    # print(x_points, y_points, hierarchy, args, kwargs)
    if len(x_points) == 0 or len(y_points) == 0:
        return x_points, y_points  # check if contour is empty
    
    if ("external_noise_size" in kwargs) and hierarchy == -1:
        external_noise_size = kwargs["external_noise_size"]
        x_points, y_points = ctr_external_denoise(x_points, y_points, noise_size=external_noise_size)
        if len(x_points) == 0 or len(y_points) == 0:
            return x_points, y_points

    if ("internal_noise_size" in kwargs) and hierarchy != -1:
        internal_noise_size = kwargs["internal_noise_size"]
        x_points, y_points = ctr_internal_denoise(x_points, y_points, noise_size=internal_noise_size)
        if len(x_points) == 0 or len(y_points) == 0:
            return x_points, y_points
    
    if "low_pass_ratio" in kwargs:
        low_pass_ratio = kwargs["low_pass_ratio"]
        x_points, y_points = ctr_fft_smooth(x_points, y_points, low_pass_ratio)

    return x_points, y_points


def ctr_internal_denoise(x_points, y_points, noise_size):
    if len(x_points) < noise_size or len(y_points) < noise_size:
        x_points = []
        y_points = []
    return x_points, y_points


def ctr_external_denoise(x_points, y_points, noise_size):
    if len(x_points) < noise_size or len(y_points) < noise_size:
        x_points = []
        y_points = []
    return x_points, y_points


def ctr_fft_smooth(x_points, y_points, low_pass_ratio=10):
    contour_len = len(x_points)
    if contour_len < 3:
        return x_points, y_points  # check if contour is too short
    
    x_points, y_points = resample_2d_contour(x_points, y_points)
    x_points, y_points = fft_low_pass(x_points, y_points, filter_ratio=low_pass_ratio)
    x_points = x_points[:contour_len]
    y_points = y_points[:contour_len]
    return x_points, y_points


def resample_2d_contour(x, y, min_pt=1024):
    repeat = math.ceil(2 + (min_pt / len(x)))
    new_x = np.tile(x, repeat)
    new_y = np.tile(y, repeat)
    return new_x, new_y


def fft_low_pass(x, y, filter_ratio=8):
    """Apply low pass filter to x and y arrays by zeroing out high freq components.
    
    ** actually, it's a band pass filter **

    Parameters
    ----------
    x, y : array_like
        Input arrays.
    length : int > len(x)
        The pass frequency length in frequency domain

    Returns
    -------
    np.real(newZ), np.imag(newZ) : array_like
        Filtered arrays of x and y
    """
    if filter_ratio >= 45:
        return x, y

    # Combine x and y into a complex signal
    z = x + 1j * y
    # Compute FFT using numpy
    C = np.fft.fft(z)

    demod = np.abs(C)
    indices = np.where(demod > 10)[0]

    wave_length = len(x)
    pass_length = int(wave_length * filter_ratio * 0.01)
    mask = np.zeros_like(C, dtype=float)
    length = pass_length + 1

    # fix the bug if filtered signal's main freq only 0
    try:
        if length < indices[1]:
            length = indices[1] + 1
    except IndexError:
        print("no other freq. > 10")
        length = int(wave_length * 10 * 0.01)
    
    mask[0] = 1
    mask[1:length] = 1
    mask[-length:] = 1

    # Apply low-pass mask
    C *= mask
    
    # Compute inverse FFT using numpy
    newZ = np.fft.ifft(C)
    return np.real(newZ), np.imag(newZ)
