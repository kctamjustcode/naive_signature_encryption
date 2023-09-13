import random, bn256, time
from bn256 import order

p = order
bng1 = bn256.curve_G
#k1, bng1 = bn256.g1_random()
bng2 = bn256.twist_G

###

a = 15287491200509037535586126781845375159542932436393982858180947378806576190056
sa = 25163806614334882268805589209959277324317058736175162994007741528647995324057
b = 51171943249810390051700320504625110284666833506464208506542823221171180356642
sb = 64206539996231185253925940673677064070322736049854092957381650183145009311451

user_A = { 'id': 'A111', 'pw': '0000', 'pri_key': a, 'send_key': sa }
user_B = { 'id': 'B111', 'pw': '0000', 'pri_key': b, 'send_key': sb }
user_list = [ user_A, user_B]

def dumb_id_to_prikey(userid: str) -> int:
    result = 0
    for user in user_list:
        if user['id'] == userid:
            result = user['pri_key']
    return result

def dumb_id_to_pubkey(userid: str) -> int:
    result = bng2
    for user in user_list:
        if user['id'] == userid:
            result = user['pri_key']*result
    assert result != bng2
    return result

def dumb_id_to_sendkey(userid: str) -> int:
    result = 0
    for user in user_list:
        if user['id'] == userid:
            result = user['send_key']
    return result


### BLS Signature
### Use py_ecc default package instead
tbls=time.time()

pub_key = dumb_id_to_prikey('A111')* bng2
sigma = dumb_id_to_prikey('A111')* bn256.g1_hash_to_point('message')
tblsn=time.time()
if bn256.optimal_ate(bng2, sigma)==bn256.optimal_ate(pub_key, bn256.g1_hash_to_point('message')):
    print('trivial BLS sign verified', tblsn-tbls)

### Aggregation
tabls=time.time()
hm1, hm2, hm3 = bn256.g1_hash_to_point('message1'), bn256.g1_hash_to_point('message2'), bn256.g1_hash_to_point('message3')
sign1 = dumb_id_to_prikey('A111')*bn256.g1_hash_to_point('message1')
sign2 = dumb_id_to_prikey('A111')*bn256.g1_hash_to_point('message2')
sign3 = dumb_id_to_prikey('A111')*bn256.g1_hash_to_point('message3')
aggre_sign = sign1+sign2+sign3

rhs = bn256.optimal_ate(bng2, aggre_sign)
lhs1 = bn256.optimal_ate(dumb_id_to_pubkey('A111'), hm1)
lhs2 = bn256.optimal_ate(dumb_id_to_pubkey('A111'), hm2)
lhs3 = bn256.optimal_ate(dumb_id_to_pubkey('A111'), hm3)
lhs = lhs1*lhs2*lhs3

#print(rhs==lhs) #manual ver.

def ecp_sum(pl: list):
    pt_sum = bn256.curve_point(bn256.gfp_1(0), bn256.gfp_1(0), bn256.gfp_1(0))
    for pt in pl:
        pt_sum = bn256.curve_point.add(pt_sum, pt)
    return pt_sum

def gfp_12_prod(vl: list):
    vf = bn256.gfp_12_one
    for item in vl:
        vf = vf* item
    return vf

hm_list = [hm1, hm2, hm3]
sign_list = [sign1, sign2, sign3]
vf_list = [lhs1, lhs2, lhs3]

rhs_f = bn256.optimal_ate(bng2, ecp_sum(sign_list))
lhs_f = gfp_12_prod(vf_list)

tablsn=time.time()
if rhs_f==lhs_f:
    print('BLS aggregated sign verified: ', tablsn-tabls)

#rhs_f_inv = bn256.gfp_12.inverse(rhs_f) #multiplicative inverse
