import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0, 3, 500)
y = x - (3/2) * x**2 + 0.5 * x**3 + 0.5*x

# Plot
plt.figure(figsize=(8,6))
plt.plot(x, y)
plt.xlabel(r'$-\frac{w}{z}$')
plt.ylabel(r'$-\frac{W}{EA}\left(\frac{L}{z}\right)^3$')
plt.title('Problem 2: Normalized force vs (-w/z)')
plt.legend()
plt.grid(True)
plt.show()
