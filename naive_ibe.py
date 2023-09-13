### IBE Encryption 3.3.1
import bn256, random, time
from bn256 import order

p = order
bng1 = bn256.curve_G
bng2 = bn256.twist_G

M = random.randint(1, p-1)
s = random.randint(1, p-1) # user A's secret_key/ hash of it
r = random.randint(1, p-1) # sender's sendkey/ hash of it

### encryption and sending from user B
def ibe_encrypt(Mesg: int, id_no: int, A_prik: int, B_sdk: int):
    M_msg = bin(Mesg)[2:] # maximal to be 256 binary bits
    M_len = len(M_msg)

    Q_id = bn256.g1_hash_to_point(str(id_no)) # public id to point; expected to be int
    U_b = B_sdk*bng2
    TxPair = bn256.optimal_ate(A_prik*bng2, B_sdk*Q_id)

    H2TPr = bn256.gt_hash(TxPair) # gt_hash = 512 binary bits
    M_H2TPr = bin(int(H2TPr.hex(),16))[2:M_len]
    V_b = bin(int(M_msg,2)^int(M_H2TPr,2))[2:]
    return (U_b, V_b, M_len)

### decryption and received by user A
def ibe_decrypt(C, id_no: int, A_prik: int):
    Q_id = bn256.g1_hash_to_point(str(id_no))
    RxPair = bn256.optimal_ate(C[0], A_prik*Q_id) # user A's secret_key

    H2RPr = bn256.gt_hash(RxPair)
    M_H2RPr = bin(int(H2RPr.hex(),16))[2: C[2]]
    return int(C[1],2)^int(M_H2RPr,2)

t331=time.time()
C = ibe_encrypt(M, 321, s, r)
D = ibe_decrypt(C, 321, s)
t331n=time.time()

if D == M:
    print("BF-IBE: ", t331n-t331)


### IBEncryption 3.5
tt=time.time()

gfp_12_one = bn256.gfp_12_one
gfp_6_zero = bn256.gfp_6_zero

alpha = random.randint(1, p-1)
g = random.randint(1, p-1)
g1 = (g*alpha)%p
g2 = random.randint(1, p-1)
g_ms = (g2*alpha)%p

U = [random.randint(1, p-1) for _ in range(256)]

#print('phrase 1 done')

v = random.randint(1, p-1)
v_bin = bin(v)[2:]
v_bin_ind = [i for i in range(len(v_bin)) if v_bin[i]=='1']
r = random.randint(1, p-1)
u_i = [U[i] for i in v_bin_ind]
sum_ui = sum(u_i)
dv = ( (g_ms+r*sum_ui)%p*bng1, (r*g)%p*bng1)

#print('phrase 2 done')

t = random.randint(1, p-1)
M = random.randint(1, 24)
M_gfp_6 = bn256.gfp_6(bn256.gfp_2(0,0), bn256.gfp_2(0,0), bn256.gfp_2(0,M)) # at least can split into 12ps
M_gfp_12 = bn256.gfp_12(gfp_6_zero, M_gfp_6)

C1 = bn256.gfp_12.exp(bn256.optimal_ate(g2*bng2, g1*bng1), t)*M_gfp_12
C2_coef = (t*g)%p
C3_coef = (t*sum_ui)%p
C = (C1, C2_coef*bng2, C3_coef*bng2)

numtr = bn256.optimal_ate(C[2], dv[1])
dentr = bn256.gfp_12.inverse(bn256.optimal_ate(C[1], (dv[0])))
DC = numtr*dentr*C[0]
#print('encrypt-decrypt sucess: ', DC==M_gfp_12)

ttn=time.time()
if M_gfp_12==DC:
    print("Wat's IBE: ", ttn-tt)
    
