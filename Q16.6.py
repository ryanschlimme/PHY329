import numpy as np
import matplotlib.pyplot as plt

def Fourier(x, i):
    a0 = 0
    f = np.array([(4/(2*n-1)/np.pi)
                   *np.sin(2*np.pi *(2*n-1)/0.25 *x) for n in range(1, i+1)])
    f = np.append(f, a0)
    return f.sum()

def FourierSingle(x, i):
    a0 = 0
    f = np.array([(4/(2*n-1)/np.pi)
                   *np.sin(2*np.pi *(2*n-1)/0.25 *x) for n in range(i, i+1)])
    f = np.append(f, a0)
    return f.sum()


xvals = np.linspace(0, 1, 1000)

y1 = np.array([FourierSingle(x,1) for x in xvals])
y2 = np.array([FourierSingle(x,2) for x in xvals])
y3 = np.array([FourierSingle(x,3) for x in xvals])
y4 = np.array([FourierSingle(x,4) for x in xvals])
y5 = np.array([FourierSingle(x,5) for x in xvals])
y6 = np.array([FourierSingle(x,6) for x in xvals])
total = np.array([Fourier(x,6) for x in xvals])

plt.plot(xvals, y1, 'r:')
plt.plot(xvals, y2, 'r:')
plt.plot(xvals, y3, 'r:')
plt.plot(xvals, y4, 'r:')
plt.plot(xvals, y5, 'r:')
plt.plot(xvals, y6, 'r:')
plt.plot(xvals, total, linewidth = 2, color = "black")
#plt.show()
plt.savefig("Q16.6.png")