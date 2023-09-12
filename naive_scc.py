import sympy, random
from sympy.abc import x, y, z

# generating indice-set, sorting indice and defining polynomials
def gen_multi_ranindlst(lngth: int, rnge: int, m: int):
    temp = set()
    while len(temp) < lngth:
        rand_elem = tuple([random.randint(0, rnge-1) for _ in range(m)])
        temp.add(rand_elem)
    temlst = list(temp)
    return temlst

def sort_on_x(tdata: list):
    xl = [tdata[i][0] for i in range(len(tdata))]
    xl.sort()
    temp = []
    for xterm in xl:
        for td in tdata:
            if td[0] == xterm:
                temp.append(td)
                tdata.remove(td)
            continue
    return temp

def double_sort(sxdata: list): #assume sorted x-coord, two coordinates
    xl = list(set([sd[0] for sd in sxdata]))
    xl.sort()
    final = []
    for xi in xl:
        tempx = []
        for sd in sxdata:
            if sd[0] == xi:
                tempx.append((sd[1],sd[0]))
        stempx = sort_on_x(tempx)
        for term in stempx:
            final.append((term[1], term[0]))
    return final

def unipoly_eval(cfl: list, pnt: int):
    m = len(cfl)
    evalsum = 0
    for i in range(m):
        evalsum += cfl[i] * pnt**[m-1 - i]
    return evalsum

def bivpoly(cfl: list, indlst: list):
    assert len(cfl) == len(indlst)
    lgth = len(cfl)
    polysum = 0
    for i in range(lgth):
        polysum += cfl[i] * x**indlst[i][0] * y**indlst[i][1]
    return polysum

def bivpoly_eval(cfl: list, indlst: list, pnt: tuple):
    assert len(cfl) == len(indlst)
    lgth = len(cfl)
    evalsum = 0
    for i in range(lgth):
        evalsum += cfl[i] * pnt[0]**indlst[i][0] * pnt[1]**indlst[i][1]
    return evalsum

def natural_sort(odata: list):
    final = []
    lgth = len(odata)
    for _ in range(lgth):
        final.append(min(odata))
        odata.remove(min(odata))
    return final

def just_triple_sort(rdata: list):
    #m = len(rdata[0][0])
    lyrs = []
    for j in range(3):
        lyr = list(set([rdata[i][j] for i in range(len(rdata))]))
        lyr.sort()
        lyrs.append(lyr)
    final = []
    for i in lyrs[0]:
        for j in lyrs[1]:
            for k in lyrs[2]:
                if (i,j,k) in rdata:
                    final.append((i,j,k))
    return final

def tripoly(cfl: list, indlst: list):
    assert len(cfl) == len(indlst)
    lgth = len(cfl)
    polysum = 0
    for i in range(lgth):
        polysum += cfl[i] * x**indlst[i][0] * y**indlst[i][1] * z**indlst[i][2]
    return polysum

def tripoly_eval(cfl: list, indlst: list, pnt: tuple):
    assert len(cfl) == len(indlst)
    lgth = len(cfl)
    evalsum = 0
    for i in range(lgth):
        evalsum += cfl[i] * pnt[0]**indlst[i][0] * pnt[1]**indlst[i][1] * pnt[2]**indlst[i][2]
    return evalsum

### Section 3
alpha = (2, 1)

coeflst = [random.randint(0, 2**13) for _ in range(2**2)]
stedind = double_sort(gen_multi_ranindlst(2**2, 2**4, 2))
ipoly = bivpoly(coeflst, stedind)

ipolyex = ipoly.as_poly().eval(alpha[0])
ipoly_eval = ipolyex.eval(alpha[1])
q1 = sympy.quo(ipoly, x - alpha[0])
r1 = sympy.rem(ipoly, x - alpha[0])
q2 = sympy.quo(r1, y - alpha[1])
r2 = sympy.rem(r1, y - alpha[1])
print('lemma 3.1 is ', r2 - ipoly_eval == 0)

beta = (3, 2, 1)
eta =  (4, 5, 6)

length = 2**3
index_range = 2**2
coeflst_tri = [random.randint(0, 2**13) for _ in range(length)]
tri = gen_multi_ranindlst(length, index_range, 3)
stri = just_triple_sort(tri)
tpoly = tripoly(coeflst_tri, stri)

tpolyex = tpoly.as_poly().eval(beta[0])
tpolyexy = tpolyex.eval(beta[1])
tpoly_eval = tpolyexy.eval(beta[2])

epolyex = tpoly.as_poly().eval(eta[0])
epolyexy = epolyex.eval(eta[1])
epoly_eval = epolyexy.eval(eta[2])

qt1 = sympy.quo(tpoly, x - beta[0])
rt1 = sympy.rem(tpoly, x - beta[0])
qt2 = sympy.quo(rt1, y - beta[1])
rt2 = sympy.rem(rt1, y - beta[1])
qt3 = sympy.quo(rt2, z - beta[2])
rt3 = sympy.rem(rt2, z - beta[2])
print('lemma 3.1 in three-variable is ', rt3 - tpoly_eval == 0)

qt1ex = qt1.as_poly().eval(eta[0])
qt1exy = qt1ex.eval(eta[1])
qt1_eval = qt1exy.eval(eta[2])
qt2exy = qt2.as_poly().eval(eta[1])
qt2_eval = qt2exy.eval(eta[2])
qt3_eval = qt3.as_poly().eval(eta[2])

import bn256
dp = bn256.order
bng1 = bn256.curve_G
bng2 = bn256.twist_G

w1 = bng1.scalar_mul(int(sympy.ZZ(qt1_eval)%dp))
w2 = bng1.scalar_mul(int(sympy.ZZ(qt2_eval)%dp))
w3 = bng1.scalar_mul(int(sympy.ZZ(qt3_eval)%dp))
FKf = bng1.scalar_mul(int(sympy.ZZ(epoly_eval)%dp))
gnv = bng1.scalar_mul(int(sympy.ZZ(-tpoly_eval)%dp))

lhs = bn256.optimal_ate(bng2, FKf+gnv)
rh1 = bn256.optimal_ate(bng2.scalar_mul(eta[0]-beta[0]), w1)
rh2 = bn256.optimal_ate(bng2.scalar_mul(eta[1]-beta[1]), w2)
rh3 = bn256.optimal_ate(bng2.scalar_mul(eta[2]-beta[2]), w3)
rhs = rh1*rh2*rh3
print('Signature of Correct Computation in 3.2 is: ', lhs == rhs)

### Section 4
def tofinitefield(q: sympy.core.numbers.Rational):
    return int(sympy.ZZ(sympy.GF(dp)(sympy.fraction(q)[0])/sympy.GF(dp)(sympy.fraction(q)[1]))%dp)

r_vec = [random.randint(1, 2**13) for _ in range(2)]

qt1_4 = sympy.quo(tpoly - tpoly_eval, r_vec[0]* (x - beta[0]) + (y - beta[1]))
rt1_4 = sympy.rem(tpoly - tpoly_eval, r_vec[0]* (x - beta[0]) + (y - beta[1]))
qt2_4 = sympy.quo(rt1_4,  r_vec[1]* (y - beta[1]) + (z - beta[2]))
rt2_4 = sympy.rem(rt1_4,  r_vec[1]* (y - beta[1]) + (z - beta[2]))
qt3_4 = sympy.quo(rt2_4,  (z - beta[2]))
rt3_4 = sympy.rem(rt2_4,  (z - beta[2]))
print('lemma 4.1 is :', 0 == rt3_4)

qt1_4ex = qt1_4.as_poly().eval(eta[0])
qt1_4exy = qt1_4ex.eval(eta[1])
qt1_4_eval = tofinitefield(qt1_4exy.eval(eta[2]))
qt2_4exy = qt2_4.as_poly().eval(eta[1])
qt2_4_eval = tofinitefield(qt2_4exy.eval(eta[2]))
qt3_4_eval = tofinitefield(qt3_4.as_poly().eval(eta[2]))

w1_4 = bng1.scalar_mul(int(sympy.ZZ(qt1_4_eval)%dp))
w2_4 = bng1.scalar_mul(int(sympy.ZZ(qt2_4_eval)%dp))
w3_4 = bng1.scalar_mul(int(sympy.ZZ(qt3_4_eval)%dp))

rh1_4 = bn256.optimal_ate(bng2.scalar_mul(r_vec[0]* (eta[0]-beta[0]) + eta[1]-beta[1]), w1_4)
rh2_4 = bn256.optimal_ate(bng2.scalar_mul(r_vec[1]* (eta[1]-beta[1]) + eta[2]-beta[2]), w2_4)
rh3_4 = bn256.optimal_ate(bng2.scalar_mul(eta[2]-beta[2]), w3_4)
print('Adaptive S.C.C. is :', lhs == rh1_4*rh2_4*rh3_4)

'''
import bn256

def convert_tuple(t: tuple, l: int): # assume pow descending from right
    m = len(t)
    l_pow = [l**(m-1-i) for i in range(m)]
    return sum([t[i]*l_pow[i] for i in range(m)])


dp = bn256.order
bng1 = bn256.curve_G
bng2 = bn256.twist_G

t1=time.time()
t_vec = tuple([random.randint(1, dp-1) for _ in range(3)])

srs = []
for i in range(3):
    for j in range(3):
        for k in range(3):
            srs.append(bng1.scalar_mul(t_vec[0]**i * t_vec[1]**j * t_vec[2]**k))


coef = [random.randint(0, dp-1) for _ in range(2**3)]
gamma = tuple([random.randint(0, 3) for i in range(3)])

triindlst = gen_multi_ranindlst(2**3, 4, 3)
til = just_triple_sort(triindlst)
bil = [(0, tile[1], tile[2]) for tile in til]
mil = [(0, 0, tile[2]) for tile in til]
tpoly = tripoly(coef, til)

#tpolyex = tpoly.as_poly().eval(gamma[0])
#tpoly_exy = tpolyex.eval(beta[1])
#tpoly_evl = tpoly_eval.eval(beta[2])
tpoly_eval = tripoly_eval(coef, til, gamma)

qt1 = sympy.quo(tpoly, x - gamma[0])
rt1 = sympy.rem(tpoly, x - gamma[0])
qt2 = sympy.quo(rt1, y - gamma[1])
rt2 = sympy.rem(rt1, y - gamma[1])
qt3 = sympy.quo(rt2, z - gamma[2])
rt3 = sympy.rem(rt2, z - gamma[2])

def ecp_sum(pl: list):
    pt_sum = bn256.curve_point(bn256.gfp_1(0), bn256.gfp_1(0), bn256.gfp_1(0))
    for pt in pl:
        assert type(pt) == bn256.curve_point
        pt_sum += pt
    return pt_sum

def scc_sign(coeflst: list, stindlst: list):
    #assert len(coeflst) == len(stindlst)
    l = len(coeflst)
    g_vec = [bng1.scalar_mul(coeflst[i]*t_vec[0]**stindlst[i][0] *t_vec[1]**stindlst[i][1] *t_vec[2]**stindlst[i][2]%dp) for i in range(l)]
    return ecp_sum(g_vec)

sig_scc = scc_sign(coef, til)
print(sig_scc)

qt1cf = qt1.as_poly().all_coeffs()[::-1]
qt1cfzmp = [int(sympy.ZZ(qc)%dp) for qc in qt1cf]
qt2cf = qt2.as_poly().all_coeffs()[::-1]
qt2cfzmp = [int(sympy.ZZ(qc)%dp) for qc in qt2cf]
qt3cf = qt3.as_poly().all_coeffs()[::-1]
qt3cfzmp = [int(sympy.ZZ(qc)%dp) for qc in qt3cf]
print(til)
print(qt1cfzmp)
print(qt2cfzmp)
print(qt3cfzmp)

sig_q1 = scc_sign(qt1cfzmp, til)
sig_q2 = scc_sign(qt2cfzmp, bil)
sig_q3 = scc_sign(qt3cfzmp, mil)
print(sig_q1)
print(sig_q2)
print(sig_q3)
'''
