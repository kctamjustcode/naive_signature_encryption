import bn256, random

p = bn256.order
bng1 = bn256.curve_G
bng2 = bn256.twist_G


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


# 4.1 The Full Signature Scheme
x = random.randint(0, p-1)
y = random.randint(0, p-1)

u = x*bng2
v = y*bng2

z = bn256.optimal_ate(bng2, bng1)

m = random.randint(1, 92)
excp = (x+m)*fp_inv(-y, p) %p

r = excp # weakly version: r == 0
while r == excp:
    r = random.randint(0, p-1)

print(r)

sigma = fp_inv(x+m+y*r, p)*bng1
print(sigma)

lhs = bn256.optimal_ate(u+m*bng2+r*v, sigma)
print(lhs==z)

