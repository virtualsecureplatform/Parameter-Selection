#!/bin/python3

import  numpy as np

# Original TFHE's parameter.

DEF_n = 500;
DEF_α = 2.44e-5;
DEF_Nbit = 10;
DEF_N = 2**DEF_Nbit;
DEF_l = 2;
DEF_Bgbit = 10;
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

print("Gate Bootstrapping Noise")
print(DEF_n*2*DEF_l*DEF_N*(DEF_β**2)*(DEF_αbk**2)+DEF_n*(1+DEF_N)*(DEF_ε**2)+DEF_N*(2**(-2*(DEF_basebit*DEF_t+1)))+DEF_t*DEF_N*(DEF_αks**2))