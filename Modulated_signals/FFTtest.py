import numpy as np
from matplotlib import pyplot as plt

samplingrate = 100
x = np.arange(0,1000,1/samplingrate)

y = np.zeros(len(x))
for i in [0.124,0.5,np.pi, 6]:
    y += i*np.sin(i*x)



fftfreq = np.fft.rfftfreq(len(y), 1/samplingrate)
fft = np.fft.rfft(y)/len(fftfreq)


plt.plot(fftfreq, np.abs(fft))
plt.axis([0,1,None,None])
plt.show()
plt.close()
