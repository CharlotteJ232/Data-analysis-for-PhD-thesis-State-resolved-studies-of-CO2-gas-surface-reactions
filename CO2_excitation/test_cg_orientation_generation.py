import numpy as np
from joblib import Parallel, delayed
import random
from matplotlib import pyplot as plt

j_range = [0,2,4,6,8,10]

folderstart = 'C:/Users/jansenc3/Surfdrive/'
folderstart = 'C:/Users/Werk/surfdrive/'
savefolder = folderstart+'Proefschrift/Rotational state distribution/Figures/241101/'
save=True

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



#http://www.cognitive-antics.net/uniform-random-orientation/
def random_projection_simple():
    return np.random.random()


#https://en.wikipedia.org/wiki/Inverse_transform_sampling
def random_theta():
    old = np.arccos(np.random.random())
    new = np.arcsin(np.random.random())
    return old

def random_rotation():
    """
    first calculates theta, the angle between the random rotation axis and the polarization axis. 
    This is done via the arccos because the distribution of theta is NOT uniform, while the distribution of the projection (the cos) is.
    """ 
    theta = random_theta()
    avg = 2/np.pi * np.sin(theta)
    return avg



#random theta
dist = Parallel(n_jobs=2)(delayed(random_theta)() for n in range(10000))
plt.title('Random theta')
plt.hist(dist,bins=10)
plt.xlabel('Theta')
plt.ylabel('Occurrence')                       
plt.show()
plt.close()
print(np.average(dist))


#circles
dist = Parallel(n_jobs=2)(delayed(random_rotation)() for n in range(10000))
plt.title('Projection of randomly oriented circles')
plt.hist(dist,bins=10)
plt.xlabel('Projection of $\mu$ on polarization axis')
plt.ylabel('Occurrence')
if save:
    plt.savefig(savefolder+'rotation.png', dpi=500)
    plt.savefig(savefolder+'rotation.pdf', dpi=500)
plt.show()
plt.close()
print(np.average(dist))

dist = Parallel(n_jobs=2)(delayed(random_projection_simple)() for n in range(10000))
plt.title('Projection of random molecule orientations')
plt.hist(dist,bins=10)
plt.xlabel('Projection of $\mu$ on polarization axis')
plt.ylabel('Occurrence')
if save:
    plt.savefig(savefolder+'norotation.png', dpi=500)
    plt.savefig(savefolder+'norotation.pdf', dpi=500)
plt.show()
plt.close()
print(np.average(dist))


#CG

    # returns cg for random m
    # < j1 m 1 0 | j1+1 m >
    # https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients#Special_cases
tableformat = 'c|'
printstring = '  '
for m in np.arange(0, 10+0.1, 1):
    printstring += '&m=$\pm$'+str(int(m))
    tableformat += 'c|'
print (printstring)
print(tableformat)
for j in [0,2,4,6,8,10]:
    printstring='j='+str(j)
    for m in np.arange(0, j+0.1, 1):
        upper = (j-m+1)*(j+m+1)
        lower = (2*j+1)*(j+1)
        printstring += '&'+str(np.round(np.sqrt(upper/lower),2))
    print(printstring+'\\\\')
    
    

# for j in j_range:
#     print(Parallel(n_jobs=2)(delayed(calc_cg)(j) for n in range(10)))