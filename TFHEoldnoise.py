#!/bin/python3

import  numpy as np
from scipy.special import erfc

# Original TFHE's parameter.

DEF_n = 500;
DEF_α = 2.44e-5;
DEF_Nbit = 10;
DEF_N = 2**DEF_Nbit;
DEF_l = 2;
DEF_Bgbit = 9;
DEF_Bg = 2**DEF_Bgbit;
DEF_αbk = 3.73e-9;
DEF_t = 8;
DEF_basebit = 2;
DEF_αks = 2.44e-5;
DEF_μ = 2**29;
DEF_ε = 1/(2*(DEF_Bg**DEF_l))
DEF_β = DEF_Bg/2

DEF_nbarbit = 11;
DEF_nbar = 2**DEF_nbarbit;
DEF_lbar = 4;
DEF_Bgbitbar = 9;
DEF_Bgbar = 2**DEF_Bgbitbar;
DEF_αbklvl02 = 2.0**-44;
DEF_tbar = 10;
DEF_basebitlvl21 = 3;
DEF_αprivks = 2**-31;
DEF_μbar = 2**61;
DEF_εbar = 1/(2*(DEF_Bgbar**DEF_lbar))
DEF_βbar = DEF_Bgbar/2

print("TFHE Gate Bootstrapping Noise")

gbnoise  = DEF_n*2*DEF_l*DEF_N*(DEF_β**2)*(DEF_αbk**2)+DEF_n*(1+DEF_N)*(DEF_ε**2)+DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+DEF_t*DEF_N*(DEF_αks**2)

print(gbnoise)

cbnoise = DEF_n*2*DEF_lbar*DEF_nbar*(DEF_βbar**2)*(DEF_αbklvl02**2)+DEF_n*(1+DEF_nbar)*(DEF_εbar**2)+DEF_nbar*(2**(-2*(DEF_basebitlvl21*DEF_tbar+1)))+DEF_tbar*DEF_nbar*(DEF_αprivks**2)

print("TFHE GB error prob")
print(erfc(1/(16*np.sqrt(2*gbnoise))))

print("TFHE Circuit Bootstrapping Noise")
print(cbnoise)

print("TFHE ROM CMUX noise")
romnoise = DEF_αbk+7*(2*DEF_l*DEF_N*(DEF_β**2)*cbnoise+(DEF_N+1)*(DEF_ε**2))+DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+DEF_t*DEF_N*(DEF_αks**2)
print(romnoise)

print("TFHE ROM error prob")
print(erfc(1/(16*np.sqrt(2*romnoise))))

print("RAM Read Noise")
rnoise = DEF_n*2*DEF_l*DEF_N*(DEF_β**2)*(DEF_αbk**2)+DEF_n*(1+DEF_N)*(DEF_ε**2)+9*(2*DEF_l*DEF_N*(DEF_β**2)*cbnoise+(DEF_N+1)*(DEF_ε**2)) + DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+DEF_t*DEF_N*(DEF_αks**2)
print(rnoise)
print("RAM Read error prob")
print(erfc(1/(16*np.sqrt(2*rnoise))))

print("Packing Switch noise")
psnoise = DEF_n*2*DEF_l*DEF_N*(DEF_β**2)*(DEF_αbk**2)+DEF_n*(1+DEF_N)*(DEF_ε**2) + DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+8*DEF_t*DEF_N*(DEF_αbk**2)
print(psnoise)


# TFHEpp's parameter.

DEF_n = 500;
DEF_α = 2.44e-5;
DEF_Nbit = 10;
DEF_N = 2**DEF_Nbit;
DEF_l = 2;
DEF_Bgbit = 9;
DEF_Bg = 2**DEF_Bgbit;
DEF_αbk = 3.73e-9;
DEF_t = 7;
DEF_basebit = 2;
DEF_αks = 2.44e-5;
DEF_μ = 2**29;
DEF_ε = 1/(2*(DEF_Bg**DEF_l))
DEF_β = DEF_Bg/2

DEF_nbarbit = 11;
DEF_nbar = 2**DEF_nbarbit;
DEF_lbar = 4;
DEF_Bgbitbar = 9;
DEF_Bgbar = 2**DEF_Bgbitbar;
DEF_αbklvl02 = 2.0**-45;
DEF_tbar = 10;
DEF_basebitlvl21 = 3;
DEF_αprivks = 2**-31;
DEF_μbar = 2**61;
DEF_εbar = 1/(2*(DEF_Bgbar**DEF_lbar))
DEF_βbar = DEF_Bgbar/2

print("TFHEpp-10ms Gate Bootstrapping Noise")
print(3*DEF_n*2*DEF_l*DEF_N*(DEF_β**2)*(DEF_αbk**2)+DEF_n*(1+DEF_N)*(DEF_ε**2)/2+DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+DEF_t*DEF_N*(DEF_αks**2))
print("TFHE Circuit Bootstrapping Noise")
print(3*DEF_n*2*DEF_lbar*DEF_nbar*(DEF_βbar**2)*(DEF_αbklvl02**2)+DEF_n*(1+DEF_nbar)*(DEF_εbar**2)/2+DEF_nbar*(2**(-2*(DEF_basebitlvl21*DEF_tbar+1)))+DEF_tbar*DEF_nbar*(DEF_αprivks**2))