import numpy as np
import matplotlib.pyplot as plt


def spectra(A=10e-3, tau=1e-6, fmax=1000, N=4001):
    """
    Compute and plot amplitude and phase spectra of x(t) = A * exp(-|t|/tau).

    Parameters:
        A   : amplitude
        tau : time constant
        fmax: maximum frequency (Hz) for plotting
        N   : number of frequency points
    """
    # Frequency axis (Hz and rad/s)
    f = np.linspace(-fmax, fmax, N)
    omega = 2 * np.pi * f

    # Analytical Fourier transform
    X = 2 * A * tau / (1 + (omega * tau) ** 2)

    # Amplitude and phase spectra
    amplitude = np.abs(X)
    phase = np.angle(X)

    # --- Plots ---
    plt.figure(figsize=(10, 6))

    # Amplitude spectrum
    plt.subplot(2, 1, 1)
    plt.plot(f, amplitude)
    plt.title(r"Amplitude Spectrum $|X(f)|$")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Amplitude")
    plt.grid(True)

    # Phase spectrum
    plt.subplot(2, 1, 2)
    plt.plot(f, phase)
    plt.title(r"Phase Spectrum $\angle X(f)$")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase [rad]")
    plt.grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    spectra(A=1.0, tau=0.2, fmax=50)
