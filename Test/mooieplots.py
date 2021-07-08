# -*- coding: utf-8 -*-
"""
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html
https://stackoverflow.com/questions/26106552/matplotlib-style-library-not-updating-when-mplstyle-files-added-deleted
https://matplotlib.org/3.2.1/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files
@author: Charlotte
"""

import numpy as np
from matplotlib import pyplot as plt

plt.style.reload_library()
plt.style.use('voorbeeld')

plt.plot(range(10), range(10))
plt.show()
plt.close()