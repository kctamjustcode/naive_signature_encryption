import hashlib, random

from py_ecc.bls12_381 import (
    add as bls12381_add,
    curve_order as bls12381_order,
    G1 as bls12381_G1,
    multiply as bls12381_multiply
)

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
    return lm % p


q = bls12381_order
P = bls12381_G1

### ECDSA
d = random.randint(1, q-1)
Q = bls12381_multiply(P, d)

m = random.randint(1, 92)

k = random.randint(1, q-1)
kP = bls12381_multiply(P, k)

r = kP[0].n
k_inv = fp_inv(k, q)
e = int(hashlib.sha256(str(m).encode()).hexdigest(), 16)
s = k_inv* (e + r*d) %q
ecdsa_sign = (r, s)

sp = fp_inv(ecdsa_sign[1], q)
evrf = int(hashlib.sha256(str(m).encode()).hexdigest(), 16)
t1 = sp* evrf %q
t2 = sp* ecdsa_sign[0] %q
R = bls12381_add(bls12381_multiply(P, t1), bls12381_multiply(Q, t2))
print('ECDSA in bls12_381 is: ', R[0].n == r)


### ECGDSA
Qp = bls12381_multiply(P, fp_inv(d, q))

s_g = d* (k*r - e) %q
ecgdsa_sign = (r, s_g)

rp = fp_inv(ecgdsa_sign[0], q)
t1_g = rp* evrf %q
t2_g = rp* ecgdsa_sign[1] %q
R_g = bls12381_add(bls12381_multiply(P, t1_g), bls12381_multiply(Qp, t2_g))
print('ECGDSA in bls12_381 is: ', R_g[0].n == r)

### ECKCDSA
hcert = random.randint(1, 2)
rk = int(hashlib.sha256(str(r).encode()).hexdigest(), 16)
m_hcert = m* hcert
ek = int(hashlib.sha256(str(m_hcert).encode()).hexdigest(), 16)
w = rk + ek %q
s_k = d* (k-w) %q
eckcdsa_sign = (rk, s_k)

wvrf = eckcdsa_sign[0] + ek
R_k = bls12381_add(bls12381_multiply(Qp, eckcdsa_sign[1]), bls12381_multiply(P, wvrf))
print('Naive ECKCDSA in bls12_381 is: ',  int(hashlib.sha256(str(R_k[0].n).encode()).hexdigest(), 16) == rk)

### bn128 ver.
from py_ecc.bn128 import (
    add as bn128_add,
    curve_order as bn128_order,
    G1 as bn128_G1,
    multiply as bn128_multiply
)

q = bn128_order
P = bn128_G1

### ECDSA
d = random.randint(1, q-1)
Q = bn128_multiply(P, d)

m = random.randint(1, 92)

k = random.randint(1, q-1)
kP = bn128_multiply(P, k)

r = kP[0].n
k_inv = fp_inv(k, q)
e = int(hashlib.sha256(str(m).encode()).hexdigest(), 16)
s = k_inv* (e + r*d) %q
ecdsa_sign = (r, s)

sp = fp_inv(ecdsa_sign[1], q)
evrf = int(hashlib.sha256(str(m).encode()).hexdigest(), 16)
t1 = sp* evrf %q
t2 = sp* ecdsa_sign[0] %q
R = bn128_add(bn128_multiply(P, t1), bn128_multiply(Q, t2))
print('ECDSA in bn128 is: ', R[0].n == r)


### ECGDSA
Qp = bn128_multiply(P, fp_inv(d, q))

s_g = d* (k*r - e) %q
ecgdsa_sign = (r, s_g)

rp = fp_inv(ecgdsa_sign[0], q)
t1_g = rp* evrf %q
t2_g = rp* ecgdsa_sign[1] %q
R_g = bn128_add(bn128_multiply(P, t1_g), bn128_multiply(Qp, t2_g))
print('ECGDSA in bn128 is: ', R_g[0].n == r)

### ECKCDSA
hcert = random.randint(1, 2)
rk = int(hashlib.sha256(str(r).encode()).hexdigest(), 16)
m_hcert = m* hcert
ek = int(hashlib.sha256(str(m_hcert).encode()).hexdigest(), 16)
w = rk + ek %q
s_k = d* (k-w) %q
eckcdsa_sign = (rk, s_k)

wvrf = eckcdsa_sign[0] + ek
R_k = bn128_add(bn128_multiply(Qp, eckcdsa_sign[1]), bn128_multiply(P, wvrf))
print('Naive ECKCDSA in bn128 is: ',  int(hashlib.sha256(str(R_k[0].n).encode()).hexdigest(), 16) == rk)
