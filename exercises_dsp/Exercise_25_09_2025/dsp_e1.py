import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = np.pi
N = 20  # Number of Fourier terms
x = np.linspace(-L, L, 1000)


# --- Define signals ---
def square_wave(x):
    return np.where(x < 0, -1, 1)


def sawtooth_wave(x):
    return x / L  # linear ramp from -1 to 1


# --- Fourier series coefficients ---
def fourier_coeffs(f, N, L):
    """Return a0, an, bn for N terms."""
    a0 = (1 / L) * np.trapezoid(f(x), x)
    a = []
    b = []
    for n in range(1, N + 1):
        an = (1 / L) * np.trapezoid(f(x) * np.cos(n * np.pi * x / L), x)
        bn = (1 / L) * np.trapezoid(f(x) * np.sin(n * np.pi * x / L), x)
        a.append(an)
        b.append(bn)
    return a0, a, b


# Compute coefficients
sq_a0, sq_a, sq_b = fourier_coeffs(square_wave, N, L)
saw_a0, saw_a, saw_b = fourier_coeffs(sawtooth_wave, N, L)

# Compute amplitudes
sq_amplitudes = [sq_a0 / 2] + [np.sqrt(an**2 + bn**2) for an, bn in zip(sq_a, sq_b)]
saw_amplitudes = [saw_a0 / 2] + [np.sqrt(an**2 + bn**2) for an, bn in zip(saw_a, saw_b)]
freqs = np.arange(0, N + 1)

# --- Plot ---
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

axs[0].stem(freqs, sq_amplitudes, basefmt=" ")
axs[0].set_title("Amplitude Spectrum: Square Wave")
axs[0].set_xlabel("Harmonic n")
axs[0].set_ylabel("Amplitude")

axs[1].stem(freqs, saw_amplitudes, basefmt=" ")
axs[1].set_title("Amplitude Spectrum: Sawtooth Wave")
axs[1].set_xlabel("Harmonic n")
axs[1].set_ylabel("Amplitude")

plt.tight_layout()
plt.show()
