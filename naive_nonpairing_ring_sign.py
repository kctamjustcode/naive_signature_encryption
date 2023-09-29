import random

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


p = 53
q = 47
n = p*q
phi_n = (p-1)*(q-1)

e = [3, 5, 7, 11] #pub keys in RSA setting
e_inv = [fp_inv(ei, phi_n) for ei in e] #skeys in RSA setting

m = random.randint(1, 92)
k = m**4 %n #simplest C:= k* sum(y) == v

v = random.randint(1, n-1)
s = random.randint(0, 3)

x = [random.randint(1, n-1) for i in range(len(e))]
y = [x[i]**e[i] %n for i in range(len(e))]
y[s] = 0

y_s =  v*fp_inv(k, n) - sum(y) %n #simplest
x[s] = y_s ** e_inv[s] %n

ring_sign = (m, e, v, x)

y_v = [x[i]**e[i] %n for i in range(len(e))]
print('naive non-pairing ring signature: ', k*sum(y_v) %n == v)
