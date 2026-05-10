import numpy as np
import matplotlib.pyplot as plt

# 1. Setup Parameters
nx = 40              # Grid points in X
nz = 80              # Vertical levels
m2_period = 44712.0  # M2 period in seconds (~12.42 hours)
n_time = 24          # Number of records to save (2 per hour)
amplitude = 0.1      # 0.1 m/s velocity

# 2. Time and Signal Generation
dt_forcing = m2_period / n_time
t = np.linspace(0, m2_period, n_time, endpoint=False)
t_hours = t / 3600.0

# Components: u (along-boundary) and v (normal)
# Using Cos/Sin creates a 90-degree phase shift (circular ellipse)
u_signal = amplitude * np.cos(2 * np.pi * t / m2_period)
v_signal = amplitude * np.sin(2 * np.pi * t / m2_period)

# 3. Create and Save Binary Files
# MITgcm expects (Time, Nz, Nx) for Southern boundary files
u_data = np.zeros((n_time, nz, nx), dtype='>f4') # Big-endian 32-bit float
v_data = np.zeros((n_time, nz, nx), dtype='>f4')

for i in range(n_time):
    u_data[i, :, :] = u_signal[i]
    v_data[i, :, :] = v_signal[i]

u_data.tofile('OBSu_m2.bin')
v_data.tofile('OBSv_m2.bin')

# 4. Plot the Phases
plt.figure(figsize=(10, 5))
plt.plot(t_hours, u_signal, 'b-o', label='u (Along-boundary / Zonal)')
plt.plot(t_hours, v_signal, 'r-s', label='v (Normal / Meridional)')
plt.axhline(0, color='black', lw=1, ls='--')
plt.xlabel('Time (hours)')
plt.ylabel('Velocity (m/s)')
plt.title(f'M2 Tidal Components: Southern Boundary (nx={nx}, nz={nz})')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

print(f"Files saved successfully.")
print(f"Set externForcingPeriod = {dt_forcing} in your 'data' file.")
print(f"Set externForcingCycle = {m2_period} in your 'data' file.")
