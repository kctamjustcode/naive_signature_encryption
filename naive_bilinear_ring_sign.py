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

def ecp_sum(pl: list):
    pt_sum = bn256.curve_point(bn256.gfp_1(0), bn256.gfp_1(0), bn256.gfp_1(0))
    for pt in pl:
        assert type(pt)==bn256.curve_point
        pt_sum += pt
    return pt_sum

def gfp_12_prod(vl: list):
    vf = bn256.gfp_12_one
    for item in vl:
        vf = vf* item
    return vf


length = 3

pri_key_list = [random.randint(1, bn256.order-1) for i in range(length)]
pub_key_list = [pri_key_list[i]* bn256.twist_G for i in range(length)]
alpha_list = [random.randint(1, bn256.order-1) for i in range(length)]

s = random.randint(0, length-1)
m = random.randint(1, 92)

sigma_s = fp_inv(pri_key_list[s], bn256.order)* (bn256.g1_hash_to_point(str(m)) \
    + (bn256.order-1)* ecp_sum([alpha_list[i]* pri_key_list[i]* bn256.curve_G for i in list(range(s))+list(range(s+1, length))]))

sigma_list = [alpha_list[i]* bn256.curve_G if i != s else sigma_s for i in range(length)]

lhs = bn256.optimal_ate(bn256.twist_G, bn256.g1_hash_to_point(str(m)))
rhs_list = [bn256.optimal_ate(pub_key_list[i], sigma_list[i]) for i in range(length)]
print('Bilinear Ring Signature is: ', lhs == gfp_12_prod(rhs_list))
