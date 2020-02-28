#!/bin/python3

import  numpy as np
from scipy.special import erfc

n = 630;
Nbit = 10;
N = 2**Nbit;
bk_l = 3;
bk_Bgbit = 7;
Bg = 2**bk_Bgbit;
bk_stdev = 2**-25;
ks_length = 8;
ks_basebit = 2;
ks_stdev = 2**-15;
μ = 2**29;
ε = 1/(2*(Bg**bk_l))
β = Bg/2

print("TFHE Gate Bootstrapping Noise Variance")

gbnoise  = n*2*bk_l*N*(β**2)*(bk_stdev**2)+n*(1+N)*(ε**2)+N*(2**(-2*(ks_basebit*ks_length+1)))+ks_length*N*(ks_stdev**2)

print(gbnoise)

print("TFHE GB error prob")
print(erfc(1/(16*np.sqrt(2*gbnoise))))

bk_Bgbit = 6;
Bg = 2**bk_Bgbit;
ε = 1/(2*(Bg**bk_l))
β = Bg/2

print("TFHE Gate Bootstrapping Noise Variance when bk_Bgbit = 6")

gbnoise  = n*2*bk_l*N*(β**2)*(bk_stdev**2)+n*(1+N)*(ε**2)+N*(2**(-2*(ks_basebit*ks_length+1)))+ks_length*N*(ks_stdev**2)

print(gbnoise)

print("TFHE GB error prob when bk_Bgbit = 6")
print(erfc(1/(16*np.sqrt(2*gbnoise))))