import math, random, sympy
from sympy.ntheory.residue_ntheory import quadratic_residues, sqrt_mod, is_quad_residue
from sympy.ntheory.modular import crt

def fp_inv(a: int, p: int) -> int:
    # Extended euclidean algorithm to find modular inverses for integers
    a %= p
    if a == 0:
        return 0
    lm, hm = 1, 0
    low, high = a % p, p
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    return lm %p

def fp_div(x, y, p:int):
    return x * fp_inv(y, p) % p

def fp_div_polys(a, b, p: int):
    """
    Long polynomial difivion for two polynomials in coefficient form
    """
    a = [x for x in a]
    o = []
    apos = len(a) - 1
    bpos = len(b) - 1
    diff = apos - bpos
    while diff >= 0:
        quot = fp_div(a[apos], b[bpos], p)
        o.insert(0, quot)
        for i in range(bpos, -1, -1):
            a[diff + i] -= b[i] * quot
        apos -= 1
        diff -= 1
    return [x % p for x in o]

def gcdExtended(a, b):
    # Base Case
    if a == 0 :
        return b,0,1

    gcd,x1,y1 = gcdExtended(b%a, a)

    # Update x and y using results of recursive
    # call
    x = y1 - (b//a) * x1
    y = x1

    return gcd, x, y


# 11.25 Rabin's Signature on HAC
def simple_rabin_sign(m, p, q):
    n = p*q
    r = random.randint(1, p*q-1)

    # simple hash by r*m
    while not is_quad_residue(r*m %n, n):
        r = random.randint(1, p*q-1)

    M = r*m %n # simplest hash
    return sqrt_mod(M, n)

def simple_rabin_verify(s, n): #problematic
    m_b = s**2 %n
    return 1 if is_quad_residue(m_b, n) else 0

# Normal Rabin's Signature
def rabin_sign(m, p , q, b):
    d = fp_div(b, 2, n)
    dp = fp_div(b, 2, p)
    dq = fp_div(b, 2, q)

    r = random.randint(1, p*q-1)
    while not is_quad_residue(r*m+d**2 %n, n):
        r = random.randint(1, p*q-1)

    c = r*m %n
    xp_p = -dp + sqrt_mod(c+dp**2, p) %n
    #xp_m = -dp - sqrt_mod(c+dp**2, p) %n
    xq_p = -dq + sqrt_mod(c+dq**2, q) %n
    #xq_m = -dq - sqrt_mod(c+dq**2, q) %n
    return int(crt([p, q], [xp_p, xq_p])[0]), r

def rabin_verify(sign, m, n, b): # simplified
    s = sign[0]
    c = s*(s+b) %n
    d = fp_div(b, 2, n)
    return c %n == sign[1]*m %n
    #return 1 if is_quad_residue(c+d**2, n) else 0


p = 53
q = 47
n = p*q

m = random.randint(1, 92)
s = simple_rabin_sign(m, p, q)
print('simple Rabin signature: ', 1 == simple_rabin_verify(s, n))

b = random.randint(1, n-1)
rs = rabin_sign(m, p, q, b)
print('Rabin Signautre: ', 1 == rabin_verify(rs, m, n, b))
