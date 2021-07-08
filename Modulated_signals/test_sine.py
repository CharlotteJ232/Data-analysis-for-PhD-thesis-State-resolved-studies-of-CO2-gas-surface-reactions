import numpy as np
from matplotlib import pyplot as plt


from unittest import TestCase

class TryTesting(TestCase):
    def test_always_passes(self):
        self.assertTrue(True)

    def test_always_fails(self):
        # arrange
        x = 1
        y = 2

        # act
        z = x - y

        # assert
        self.assertGreater(z, 0)

def generate_sine():
    amp = 0.1
    t = np.arange(0,5,0.1)
    t2 = t + 256256
    f = 0.5

    a = np.sin(t*f*2*np.pi)
    b = np.sin(t2*f*2*np.pi)

    return (a-b)**2
