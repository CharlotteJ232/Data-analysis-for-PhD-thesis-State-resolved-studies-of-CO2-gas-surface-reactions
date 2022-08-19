import numpy as np
from matplotlib import pyplot as plt

def TOF_function_time(t, L, v0, ts, alpha):
    return (L/(t-ts))**4 * np.exp(-((L/(t-ts)-v0)/alpha)**2)


t = np.arange(0, 2, 0.0001)
ts = 0.40001
L = 760.8
alpha = 176
v0 = 1568

normal = TOF_function_time(t, L, v0, ts, alpha)
small_L = TOF_function_time(t, 1, v0, ts, alpha)
# plt.plot(t, normal, label='normal')
# plt.plot(t, small_L, label='small L')
for L in np.arange(700.5, 0, -100):
    fun = TOF_function_time(t, L, v0, ts, alpha)
    plt.scatter(t[np.argmax(fun)], np.max(fun), color='gray')
    plt.scatter(np.sum(t*fun)/np.sum(fun), np.max(fun), color='black')
    plt.plot(t, fun, label=str(L))
plt.legend()
plt.show()
plt.close()

print(t[np.argmax(small_L)])


fig, ax = plt.subplots(figsize=(15,12))
first = True
for L in np.arange(700.5, 0, -50):
    fun = TOF_function_time(t, L, v0, ts, alpha)
    if first:
        plt.plot([0,L], [ts,t[np.argmax(fun)]],color='gray')
        plt.plot([0,L], [ts,np.sum(t*fun)/np.sum(fun)], color='black')
        first=False
    plt.plot(L, t[np.argmax(fun)], 'o',color='gray')
    plt.plot(L, np.sum(t*fun)/np.sum(fun), 'o', color='black')
plt.scatter(-10, -10, color='gray', label='Most probable')
plt.scatter(-10, -10, color='black', label='Average')
plt.ylabel('Time (ms)')
plt.xlabel('Path length (mm)')
plt.title('Most probable vs average')
plt.legend()
plt.axis([0, None, ts, 0.9])
plt.show()
plt.close()

# for L in np.arange(700, 20, -50):
#     fun = TOF_function_time(t, L, v0, ts, alpha)
#     plt.scatter(L, np.sum(t*fun)/np.sum(fun))
# plt.scatter(0, ts)
# plt.ylabel('Average time')
# plt.xlabel('Path length')
# plt.title('Average')
# plt.axis([0, None, ts, None])
# plt.show()
# plt.close()

#conclusion: extrapolating maximum in time domain to L=0 will result in ts

