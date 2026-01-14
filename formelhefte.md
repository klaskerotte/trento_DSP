# Formulasheet digital signal processing

## Infinite geometric series sum

$$
\sum_{k=0}^{\infty}(ar^k) = \frac{a}{1-r}
$$

## Finite geometric series sum 

$$
\sum_{k=0}^{n} (ar^k) = \frac{ar^{n+1}-a}{r-1}
$$

## DTFT and IDFT 

$$
X(j\omega) = \sum{-\infty}^{\infty}(x[n]e^{-j\omega n}
$$

$$
x[n] = \frac{1}{2\pi}\int_{-\pi}^{\pi}X(j\omega)e^{j\omega n}~d\omega
$$

## N-point DFT and IFT

$$
x[k] = \sum_{n = 0}^{N-1} x[n]e^{-j\frac{2\pi n k}{N}}
$$

$$
x[n] = \frac{1}{2\pi} \sum_{k = 0}^{N-1}x[k]e^{j\frac{2\pi n k}{N}}
$$


## Common z-transformations

For causal (rightsided) systems:

$$
z^{-k} = \delta [n -k]
$$

$$
\sum_{k=0}^{\infty}(z^{-1})^{k} = \frac{1}{1-z^{-1}} = u[n] 
$$


$$
\frac{1}{1-z^{-1}} \cdot z^{-k} = u[n-k]
$$

For anticausal (leftsided) systems replace $(n)$ by $(-n)$

## Bilinear transform

s = \frac{2}{T}\frac{1-z^{-1}}{1+z^{-1}}

## Pole-zero shenanigans

For causal systems:
- Look outwards from the outermost pole
- If the region of convergence (ROC) STRICTLY contains the unit circle, then the system is stable.

For anticausal systems:
- Look inwards from the innermost pole
- If the ROC STRICTLY contains the unit circle, then the system is stable.

For twosided systems:
- Look between two poles.
- If the ROC STRICTLY contains the unit circle, then the system is stable.

## Filters

### Causal stable low pass filter
- Pole in $z = 0$
- Zero in $z = -1$

$$
H(z) = 1 + z^{-1} = \frac{z+1}{z} 
$$

The low pass qualities can be shown by examining the frequency response at $\omega = 0$ and $\omega = \pi$.

$$
H(j\omega) = \frac{e^{j\omega} + 1}{e^{j\omega}} \to H(j0) = 2 ~ \cup ~ H(j\pi) = 0
$$

ROC for this particular filter is $|z| > 0$, meaning the whole plane except the origin.


### Causal stable high pass filter
- Pole in $z = 0$
- Zero in $z = 1$

$$
H(z) = 1 - z^{-1} = \frac{z - 1}{z} 
$$

The high pass qualities can be shown by examining the frequency response at $\omega = 0$ and $\omega = \pi$.


$$
H(j\omega) = \frac{e^{j\omega}-1}{e^{j\omega}} \to H(j0) = 0 ~ \cup ~ H(j\pi) = 2
$$

The ROC for this filter is $|z| > 0$, meaning the whole plane except the origin.

### Causal unstable band pass filter

- The angle of the complex zeros determine the passband
- The complex zeros come as conjugate pairs in practical scenarios (nonconjugate pairs are theoretically possible)

$$
H(z) = \frac{z^2 - 1}{z^2} = \frac{(z+1)(z-1)}{z^2}
$$

$$
H(j\omega) = \frac{e^{j2\omega} - 1}{e^{j2\omega}} \to H(j0) = 0 ~ \cup ~ H(j\pi) = 0 ~ \cup ~ 
H(j\frac{\pi}{2}) = 2
$$


### Causal unstable band stop filter

- The angle of the complex zeros determine the stopband
- The complex zeros come as conjugate pairs in practical scenarios (nonconjugate pairs are theoretically possible)


$$
H(z) = \frac{z^2 + 1}{z^2} = \frac{(z+j)(z-j)}{z^2}
$$

$$
H(j\omega) = \frac{e^{j2\omega} + 1}{e^{j2\omega}} \to H(j0) = 2 ~ \cup ~ H(j\pi) = 2 ~ \cup ~ 
H(j\frac{\pi}{2}) = 0
$$

