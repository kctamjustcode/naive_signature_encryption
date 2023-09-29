import bn256, random

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


### identity-based
q = bn256.order
bng1 = bn256.curve_G
bng2 = bn256.twist_G


s = random.randint(1, q-1)
Ppub = s*bng2
id = random.randint(1, 92)
Qid = bn256.g1_hash_to_point(str(id))
Sid = s*Qid


r = random.randint(1, q-1) # randomness in multiple messages pass to alpha
U = r*Qid # signer's granted blinding

m = random.randint(1, 92)
alpha = random.randint(1, q-1) # blinding 1
beta = random.randint(1, q-1) # blinding 2

Up = alpha*U + alpha*beta*Qid
hm = m* bn256.g1_compress(Up)[0]
h = fp_inv(alpha, q)* hm + beta


V = (r+h)*Sid

Vp = alpha*V # unblinding 1

print( bn256.optimal_ate(bng2, Vp) == bn256.optimal_ate(Ppub, Up+hm*Qid) )


### Aggregating Signatures
m2 = random.randint(1, 92)
alpha2 = random.randint(1, q-1)
beta2 = random.randint(1, q-1)

Up2 = alpha2*U + alpha2*beta2*Qid
hm2 = m2* bn256.g1_compress(Up2)[0]
h2 = fp_inv(alpha2, q)* hm2 + beta2


m3 = random.randint(1, 92)
alpha3 = random.randint(1, q-1)
beta3 = random.randint(1, q-1)

Up3 = alpha3*U + alpha3*beta3*Qid
hm3 = m3* bn256.g1_compress(Up3)[0]
h3 = fp_inv(alpha3, q)* hm3 + beta3


V2 = (r+h2)*Sid
Vp2 = alpha2*V2

V3 = (r+h3)*Sid
Vp3 = alpha3*V3

#print( bn256.optimal_ate(bng2, Vp2) == bn256.optimal_ate(Ppub, Up2+hm2*Qid) )
#print( bn256.optimal_ate(bng2, Vp3) == bn256.optimal_ate(Ppub, Up3+hm3*Qid) )
print( bn256.optimal_ate(bng2, Vp+Vp2+Vp3) == bn256.optimal_ate(Ppub, (Up+Up2+Up3) + (hm+hm2+hm3)* Qid) )
