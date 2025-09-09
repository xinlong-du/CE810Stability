#CE 810 HW 2, Pragun Shrestha
import numpy as np
import matplotlib.pyplot as plt
#Problem 1

# Define x range (normalized deflection)
x = np.linspace(0, 3, 300)

# Corrected normalized load expression
y = x - (3/2)*x**2 + (1/2)*x**3

# Plot
plt.figure(figsize=(7,5))
plt.plot(x, y, 'b-', linewidth=2)

plt.xlabel(r'Normalized Deflection, $-w/z$')
plt.ylabel(r'Normalized Load, $-W(L/z)^3/[EA]$')
plt.title('Normalized Load-Deflection Relationship')

plt.xlim(0, 3)
plt.ylim(-0.4, 0.4)
plt.axhline(0, color='k', linewidth=0.8)  # x-axis
plt.axvline(0, color='k', linewidth=0.8)  # y-axis
plt.grid(True)
plt.show()