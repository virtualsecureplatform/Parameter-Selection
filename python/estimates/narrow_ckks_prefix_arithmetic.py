"""Independent integer phase arithmetic for private CKKS diagnostic prefixes.

Fifteen-bit limbs and one 31-bit scratch prime give exact signed convolutions
for N <= 16384. This does not use TFHEpp arithmetic or serialize key material.
"""
import cmath
import math

PRIME = 2013265921
GENERATOR = 31
DIGIT_BITS = 15


def ntt(values, inverse=False):
    values = [x % PRIME for x in values]
    n = len(values)
    if n < 2 or n & (n-1) or (PRIME-1) % n:
        raise ValueError('unsupported transform length')
    root = pow(GENERATOR, (PRIME-1)//n, PRIME)
    if pow(root, n, PRIME) != 1 or pow(root, n//2, PRIME) == 1:
        raise ValueError('scratch root has incorrect order')
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            values[i], values[j] = values[j], values[i]
    length = 2
    while length <= n:
        step = pow(GENERATOR, (PRIME-1)//length, PRIME)
        if inverse:
            step = pow(step, PRIME-2, PRIME)
        half = length//2
        for start in range(0, n, length):
            weight = 1
            for i in range(start, start+half):
                left, right = values[i], values[i+half]*weight % PRIME
                values[i] = (left+right) % PRIME
                values[i+half] = (left-right) % PRIME
                weight = weight*step % PRIME
        length *= 2
    if inverse:
        normalizer = pow(n, PRIME-2, PRIME)
        values = [x*normalizer % PRIME for x in values]
    return values


class PhaseDecoder:
    def __init__(self, secret):
        self.secret = list(secret)
        self.n = len(secret)
        if self.n < 2 or self.n > 16384 or self.n & (self.n-1):
            raise ValueError('unsupported ring dimension')
        if any(x not in (-1, 0, 1) for x in secret):
            raise ValueError('phase decoder needs a ternary secret')
        if self.n*((1 << DIGIT_BITS)-1) >= PRIME//2:
            raise ValueError('signed convolution exceeds scratch prime budget')
        self.spectrum = None

    def phase(self, a, b, log_q, force_ntt=False):
        if len(a) != self.n or len(b) != self.n or not 1 <= log_q <= 448:
            raise ValueError('invalid ciphertext shape or level')
        modulus = 1 << log_q
        a = [x % modulus for x in a]
        product = [0]*self.n
        if self.n <= 512 and not force_ntt:
            for j, sign in enumerate(self.secret):
                if sign:
                    for i, value in enumerate(a):
                        position = i+j
                        product[position % self.n] += sign*value*(1 if position < self.n else -1)
        else:
            if self.spectrum is None:
                self.spectrum = ntt(self.secret+[0]*self.n)
            mask = (1 << DIGIT_BITS)-1
            for offset in range(0, log_q, DIGIT_BITS):
                spectrum = ntt([(x >> offset) & mask for x in a]+[0]*self.n)
                convolution = ntt([x*y % PRIME for x, y in zip(spectrum, self.spectrum)], True)
                for i in range(self.n):
                    value = (convolution[i]-convolution[i+self.n]) % PRIME
                    if value > PRIME//2:
                        value -= PRIME
                    product[i] += value << offset
        return [(y-x+modulus//2) % modulus-modulus//2 for x, y in zip(product, b)]


def canonical_decode(phase, scale, direct=False):
    n = len(phase)
    coefficients = [math.ldexp(float(x), -scale) for x in phase]
    indices = [(pow(5, i, 2*n)-1)//2 for i in range(n//2)]
    if direct:
        slots = []
        for index in indices:
            root = cmath.exp(1j*math.pi*(2*index+1)/n)
            value = 0j
            for coefficient in reversed(coefficients):
                value = value*root+coefficient
            slots.append(value)
        return slots
    values = [x*cmath.exp(1j*math.pi*j/n) for j, x in enumerate(coefficients)]
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            values[i], values[j] = values[j], values[i]
    length = 2
    while length <= n:
        step = cmath.exp(2j*math.pi/length)
        half = length//2
        for start in range(0, n, length):
            weight = 1+0j
            for i in range(start, start+half):
                left, right = values[i], values[i+half]*weight
                values[i], values[i+half] = left+right, left-right
                weight *= step
        length *= 2
    return [values[index] for index in indices]


def canonical_component_decode(phase, scale, direct=False):
    """Separate real/imaginary phases exactly before floating-point decoding.

Conjugation sends X to X^-1; X^(n/2) evaluates to i at every canonical
slot. Keep twice each component as integers so dividing by two cannot round
away a small component next to a much larger one. Phases are centered integer
coefficients, as returned by PhaseDecoder.phase.
"""
    n = len(phase)
    if n < 4 or n & (n-1):
        raise ValueError('unsupported component phase dimension')
    real_twice = [2*phase[0]] + [phase[i]-phase[n-i] for i in range(1,n)]
    imaginary_twice = [0] + [phase[i]+phase[n-i] for i in range(1,n)]
    half = n//2
    imaginary_twice = imaginary_twice[half:] + [-x for x in imaginary_twice[:half]]
    real = canonical_decode(real_twice,scale+1,direct)
    imaginary = canonical_decode(imaginary_twice,scale+1,direct)
    return [complex(x.real,y.real) for x,y in zip(real,imaginary)]
