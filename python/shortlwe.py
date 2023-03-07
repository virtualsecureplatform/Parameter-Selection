#!/bin/python3
import  numpy as np
from scipy.special import erfc

q = 2**8
noise = q * 2**-2

print(erfc((q/4)/np.sqrt(2*(noise))))
