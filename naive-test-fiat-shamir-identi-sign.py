import math, random, sympy
from sympy.ntheory.residue_ntheory import quadratic_residues, sqrt_mod

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

def poly_eval_modn(poly: list, a: int, n: int):
    return sum([poly[i]*a**i %n for i in range(len(poly))])%n


p = 47
q = 43
n = p*q


# identification scheme
I = 'AaBb'
i_bin = bin(int(I.encode().hex(), 16))[2:]
f_I = [int(coef) for coef in i_bin[::-1]] #naive f_I
#f_IR = [random.randint(1, n-1)* int(coef) for coef in i_bin[::-1]]

v = [poly_eval_modn(f_I, j, n) for j in range(len(f_I))]
qrn = quadratic_residues(n)
qrn.remove(0) #not absolute
qrn.remove(1)

qrv = [(i, v[i]) for i in range(len(v)) if v[i] in qrn]
vinv = [fp_inv(qrvi[1], n) for qrvi in qrv]
sqrtvinv = [min(sqrt_mod(vii, n, True)) for vii in vinv]
s = sqrtvinv
s_ind = [qr[0] for qr in qrv]

#clean_zero_by_multiple_of_primes
s_czmi = [(s_ind[i], s[i], qrv[i][1]) for i in range(len(qrv)) if qrv[i][1]*(s[i])**2 %n != 0]


t = len(f_I)
k = len(s_czmi)

for i in range(t):
    r = random.randint(0, n-1)
    x = r**2 %n
    rbv = [random.randint(0, 1) for _ in range(k)]
    si = 1
    vi = 1

    for j in range(len(rbv)):
        if rbv[j] != 0:
            si *= s_czmi[j][1]
            vi *= s_czmi[j][2]
            assert s_czmi[j][2]* s_czmi[j][1]**2 %n == 1

    yi = r*si %n

    fsbool = True
    if not x == yi**2 *vi %n:
        print('bugs', r, x, yi)
        fsbool = False

print('Fiat-Shamir Identification Scheme is: ', fsbool)


# signature scheme
r_vec = [random.randint(0, n-1) for _ in range(t)]
x_vec = [ri**2 %n for ri in r_vec]

m = random.randint(1, 92)
m_bin = bin(m)[2:]
k_m = min(k, len(m_bin))
f_m = [int(m_bin[len(m_bin)-1-i]) for i in range(k_m)]

e = []
for i in range(t):
    for j in range(k_m):
        temp = bin(poly_eval_modn(f_m, x_vec[i], n))[2:]
        e.append([int(temp[len(temp)-1-i]) for i in range(len(temp))][:k_m])

y = []
for i in range(t):
    si = r_vec[i]
    for j in range(k_m):
        if e[i][j] == 1:
            si *= s_czmi[j][1]
    y.append(si %n)

z = []
for i in range(t):
    zi = y[i]**2
    for j in range(k_m):
        if e[i][j] == 1:
            zi *= s_czmi[j][2]
    z.append(zi %n)

ze = []
for i in range(t):
    for j in range(k_m):
        temp = bin(poly_eval_modn(f_m, z[i], n))[2:]
        ze.append([int(temp[len(temp)-1-i]) for i in range(len(temp))][:k_m])

print('Fiat-Shamir Signature Scheme is: ', ze == e)
# or, assert z[i] == e[i] directly
# print([z[i] == x_vec[i] for i in range(t)])
