#!/bin/python3

import  numpy as np
from scipy.special import erfc
import gmpy2
from gmpy2 import mpfr
from noiseestimation.keyvariation import *

gmpy2.get_context().precision=200

# from noiseestimation.params.λ128bit import *
from noiseestimation.params.concrete import *

print("Gate")
print("BR noise")
brnoise = brnoisecalc(lvl01param)
print(brnoise)
print(np.sqrt(brnoise)/lvl1param.q)
print(erfc((lvl1param.q/8)/np.sqrt(2*2*brnoise)))

iksnoise = iksnoisecalc(lvl10param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/lvl0param.q)
print(erfc((lvl0param.q/8)/np.sqrt(2*iksnoise)))

print("BR Round Noise")
roundnoise = brroundnoise(lvl0param,lvl1param)
print(roundnoise)
print(np.sqrt(roundnoise)/lvl0param.q)
print(erfc((lvl0param.q/8)/np.sqrt(2*roundnoise)))

print("Gate Error")
gatenoise = (roundnoise+iksnoise+2*brnoise*(lvl0param.q/lvl1param.q)**2)
print(np.sqrt(gatenoise)/lvl0param.q)
print(erfc((lvl0param.q/8)/np.sqrt(2*gatenoise)))

print("m = 2 BR noise")
brnoise = unrollbrnoisecalc(lvl0param,lvl1param,2)
print(brnoise)
print(np.sqrt(brnoise)/lvl0param.q)
print(erfc((lvl2param.q/16)/np.sqrt(2*brnoise)))

print("m=2 Gate Error")
print(erfc((lvl0param.q/16)/np.sqrt(2*(iksnoise+brnoise))))

print("lvl02 Gate")
print("BR noise")
brnoise = brnoisecalc(lvl02param)
print(brnoise)
print(np.sqrt(brnoise)/lvl2param.q)
print(erfc((lvl2param.q/16)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl20param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/lvl0param.q)
print(erfc((lvl0param.q/16)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(brnoise*(lvl0param.q/lvl2param.q)**2)
print(erfc((lvl0param.q/16)/np.sqrt(2*(iksnoise+brnoise*(lvl0param.q/lvl2param.q)**2))))

print("BRlvl02")
brnoise = brnoisecalc(lvl02param)
print(brnoise)
print(np.sqrt(brnoise)/lvl2param.q)
print(erfc((lvl2param.q/16)/np.sqrt(2*brnoise)))
brnoise *= ((lvl1param.q/lvl2param.q)**2)

privksnoise = privksnoisecalc(lvl21param)
print("PrivIKS noise")
print(privksnoise)
print(np.sqrt(privksnoise)/lvl1param.q)
print(erfc((lvl1param.q/16)/np.sqrt(2*privksnoise)))

print("CB round noise")
roundnoise = brroundnoise(lvl0param,lvl2param,ϑ=np.ceil(np.log2(lvl2param.l)))
print(roundnoise)
print(np.sqrt(roundnoise)/lvl0param.q)
print(erfc((lvl0param.q/4)/np.sqrt(2*(roundnoise+iksnoise+2*brnoise*(lvl0param.q/lvl2param.q)**2))))

print("CB")

print(np.sqrt(cbnoisecalc(lvl02param,lvl21param))/lvl1param.q)
print(erfc((lvl1param.q/16)/np.sqrt(2*(cbnoisecalc(lvl02param,lvl21param)))))

print("MRLWE IKS noise")
mrlweiks = mrlweikscalc(lvl21mrlweparam,lvl2param)
print(mrlweiks)
print(privksnoise)
print(erfc((lvl0param.q/16)/np.sqrt(2*mrlweiks)))

print(brnoise+privksnoise-cbnoisecalc(lvl02param,lvl21param))

ROMaddress = 7 # 4 word block
print("TFHE ROM CMUX noise")
romnoise = romnoisecalc(lvl02param,lvl21param,ROMaddress)
print(romnoise)
print(np.sqrt(romnoise)/lvl1param.q)
print("TFHE ROM error prob")
print(erfc(lvl1param.q/(4*np.sqrt(2*romnoise))))

brnoise = brnoisecalc(lvl01param)

privksnoise = privksnoisecalc(lvl11param)
print("PrivIKS lvl11 noise")
print(privksnoise)
print(np.sqrt(privksnoise)/lvl1param.q)
print(erfc((lvl1param.q/16)/np.sqrt(2*privksnoise)))

print("CBlvl01")
print(erfc((lvl1param.q/16)/np.sqrt(2*(privksnoise+brnoise))))
print("CBlvl01 extp")
extpnoise = extpnoisecalc(lvl1param,privksnoise+brnoise,0,1,0)#+brnoisecalc(lvl0param,lvl1param)
print(erfc((lvl1param.q/16)/np.sqrt(2*(extpnoise))))

print("Annihilate CB lvl22")
brnoise = brnoisecalc(lvl02param)
annihilatenoise = annihilatecalc(lvl2param,brnoise)
print(np.sqrt(annihilatenoise)/lvl2param.q)
annihilatecbnoise = extpnoisecalc(lvl2param,lvl2param.σ,annihilatenoise,lvl2param.expectation_key_coefficient,lvl2param.variance_key_coefficient)
print(erfc((lvl2param.q/16)/np.sqrt(2*(annihilatecbnoise))))

print("TFHE ROM CMUX noise lvl22")
romnoise = romnoisecalc(lvl02param,lvl22param,ROMaddress)
print(romnoise)
print(np.sqrt(romnoise)/lvl2param.q)
print("TFHE ROM error prob lvl22")
print(erfc(lvl2param.q/(4*np.sqrt(2*romnoise))))

print("TFHE ROM Annihilate CMUX noise lvl22")
romnoise = annihilateromnoisecalc(lvl02param,ROMaddress)
print(romnoise)
print(np.sqrt(romnoise)/lvl2param.q)
print("TFHE Annihilate ROM error prob lvl22")
print(erfc(lvl2param.q/(4*np.sqrt(2*romnoise))))

RAMaddress = 7
RAMwordbit = 8

print("RAM Read Noise")
rnoise = ramnoisecalc(lvl01param,lvl10param,lvl02param,lvl21param,RAMaddress)
print(rnoise)
print("RAM Read error prob")
print(erfc(1/(16*np.sqrt(2*rnoise))))
print("RAM Read error prob in 3000 cycle")
print(1-(1-mpfr(erfc(1/(16*np.sqrt(2*rnoise)))))**(3000*RAMwordbit))

print("multi-bit LUT")
print("BR noise")
brnoise = brnoisecalc(lvl02param)
print(brnoise)
print(np.sqrt(brnoise)/lvl2param.q)
print(erfc((lvl2param.q/64)/np.sqrt(2*brnoise)))

iksnoise = iksnoisecalc(lvl20param)
print("IKS noise")
print(iksnoise)
print(np.sqrt(iksnoise)/lvl2param.q)
print(erfc((lvl0param.q/64)/np.sqrt(2*iksnoise)))

print("Gate Error")
print(erfc((lvl0param.q/64)/np.sqrt(2*(iksnoise+brnoise))))
