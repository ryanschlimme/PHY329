import numpy as np
import matplotlib.pyplot as plt

def T(x):
    return 3.017*np.exp(np.sqrt(0.15)*x)+236.983*np.exp(-np.sqrt(0.15)*x)

x = np.linspace(0, 10, 1000)

plt.figure(figsize=(4,3))
plt.plot(x, T(x))
plt.xlim(0,10)
plt.ylim(45,260)
plt.savefig("24.1a.png")
