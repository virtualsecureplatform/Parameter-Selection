#!/bin/python3

import  numpy as np
from scipy.special import erfc
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=200

# 128bit TFHE's parameter.
lvl0_μ = lvl1_μ = 2**29;
lvl2_μ = 2**61;

class lvl0param:
    n = 635
    α = 2**-15

class lvl1param:
    nbit = 10
    n = 2**nbit
    l = 3
    Bgbit = 6
    Bg = 2**Bgbit
    α = 2**-25
    ε = 1/(2*(Bg**l))
    β = Bg/2

class Annihilatelvl1param:
    nbit = 10
    n = 2**nbit
    l = 12
    Bgbit = 2
    Bg = 2**Bgbit
    α = 2**-25
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl2param:
    nbit = 11
    n = 2**nbit
    l = 4
    Bgbit = 9
    Bg = 2**Bgbit
    α = 2**-44
    ε = 1/(2*(Bg**l))
    β = Bg/2

class lvl10param:
    t = 7
    basebit = 2

class lvl21param:
    t = 8
    basebit = 3

class lvl22param:
    t = 10
    basebit = 3

class lvl20param:
    t = 7
    basebit =  2

