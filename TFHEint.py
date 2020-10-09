#!/bin/python3

import  numpy as np
from scipy.special import erfc
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=200

# 128bit TFHE's parameter.

DEF_n = 634;
DEF_α = 2**-15;
DEF_Nbit = 10;
DEF_N = 2**DEF_Nbit;
DEF_nbarbit = 11;
DEF_nbar = 2**DEF_nbarbit;
DEF_ltilde = 3;
DEF_Bgbittilde = 11;
DEF_Bgtilde = 2**DEF_Bgbittilde;
DEF_ttilde = 7;
DEF_basebitlvl20 = 2;
DEF_αkslvl20 = DEF_α;
DEF_αbklvl02 = 2.0**-44;

DEF_εtilde = 1/(2*(DEF_Bgtilde**DEF_ltilde))
DEF_βtilde= DEF_Bgtilde/2

print("TFHE Gate Bootstrapping TLWE2TLLE lvl02 Noise")

gbl2lnoise  = DEF_n*2*DEF_ltilde*DEF_nbar*(DEF_βtilde**2)*(DEF_αbklvl02**2)+DEF_n*(1+DEF_nbar)*(DEF_εtilde**2)

print(gbl2lnoise)

print("TFHE Gate Bootstrapping Noise")

ksnoise = DEF_nbar*(2**(-2*(DEF_basebitlvl20*DEF_ttilde+1)))+DEF_ttilde*DEF_nbar*(DEF_αkslvl20**2)
gbnoise = gbl2lnoise + ksnoise
print(gbnoise)

print("9 bit int Bootstrapping error prob")
print(erfc(1/((512)*np.sqrt(2*gbl2lnoise))))

print("Add noise")
print(2*8*gbl2lnoise)

print("Bootstrapping lvl 120")
print()

print("Signigicant term")
print(DEF_n*2*DEF_ltilde*DEF_nbar*(DEF_βtilde**2)*(DEF_αbklvl02**2))
print(DEF_n*(1+DEF_nbar)*(DEF_εtilde**2))
print(DEF_nbar*(2**(-2*(DEF_basebitlvl20*DEF_ttilde+1))))
print(DEF_ttilde*DEF_nbar*(DEF_αkslvl20**2))