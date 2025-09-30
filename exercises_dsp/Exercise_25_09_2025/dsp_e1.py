import numpy as np
import matplotlib.pyplot as plt

# Parameters
T = 4.0  # period in seconds
t_max = 4.0  # total time to plot
dt = 0.001  # time step
delay = 2.0  # delay in seconds

# Time vector
t = np.arange(0, t_max, dt)


# Delayed square wave: shift time by 'delay'
y_delayed = np.where((((t + 4) - delay) - 2 % T) < (T / 2), 1, 0) * np.exp(-2 * (t + 4))

# Plot
plt.figure(figsize=(8, 4))
# plt.step(t, y, where="post", label="Original")
plt.step(t, y_delayed, where="post", label=f"Delayed by {delay}s")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("v(t) = s(t+4)")
# plt.ylim(0, 1.5)
plt.grid(True)
plt.legend()
plt.show()
