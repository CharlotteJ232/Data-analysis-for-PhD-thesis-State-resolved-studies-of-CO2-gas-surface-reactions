import numpy as np
from joblib import Parallel, delayed
import random
from matplotlib import pyplot as plt

j_range = [0,2,4,6,8,10]

#< j1 m 1 0 | j1+1 m >

def calc_cg(j):
    m_list = np.arange(-j, j+0.1, 1)
    m = random.choice(m_list)
    upper = (j-m+1)*(j+m+1)
    lower = (2*j+1)*(j+1)
    print(j)
    print(np.sqrt(upper/lower))
    return(m)

def random_projection():
    """
    This just returns a uniform distribution
    """
    x = np.random.randn()
    y = np.random.randn()
    z = np.random.randn()
    return(z / np.sqrt((x**2+y**2+z**2)))

def random_projection_simple():
    return np.random.random()


def random_rotation():
    """
    """
    theta = np.arccos(np.random.random())
    avg = 2/np.pi * np.sin(theta)
    return avg


dist = Parallel(n_jobs=2)(delayed(random_rotation)() for n in range(10000))
plt.title('Distribution of averaged circles with random orientations')
plt.hist(dist,bins=10)
plt.xlabel('Projection of mu on polarization axis')
plt.ylabel('Occurrence')
plt.show()
plt.close()
print(np.average(dist))

dist = Parallel(n_jobs=2)(delayed(random_projection_simple)() for n in range(10000))
plt.title('Distribution of random molecule orientations')
plt.hist(dist,bins=10)
plt.xlabel('Projection of mu on polarization axis')
plt.ylabel('Occurrence')
plt.show()
plt.close()
print(np.average(dist))

    
    

# for j in j_range:
#     print(Parallel(n_jobs=2)(delayed(calc_cg)(j) for n in range(10)))