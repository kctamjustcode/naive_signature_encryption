import random, bn256, time
from bn256 import order

p = order
bng1 = bn256.curve_G
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

#####
### ElGamal
telgl=time.time()

m = random.randint(1, p-1)
M = m*bng1

A_prk = dumb_id_to_prikey('A111')
B_sk = dumb_id_to_sendkey('B111')

def encrypt_elgamal(A_priky: int, B_sdk: int, message_point: tuple) -> tuple:
    base_point = bn256.g1_hash_to_point(str(A_priky))
    receiver_point = A_priky*base_point
    cipher = B_sdk*receiver_point + bn256.g1_uncompress(message_point)
    return bn256.g1_compress(cipher)

def decrypt_elgamal(A_priky: int, B_sdk: int, cipher_point: tuple) -> tuple:
    base_point = bn256.g1_hash_to_point(str(A_priky))
    sender_point = A_priky*base_point
    decipher_point = B_sdk*sender_point
    recovered_M = bn256.g1_uncompress(cipher_point) + (bn256.order-1)*decipher_point
    return bn256.g1_compress(recovered_M)

encrypted_message = encrypt_elgamal(A_prk, B_sk, bn256.g1_compress(M))
decrypted_message = decrypt_elgamal(A_prk, B_sk, encrypted_message)

telgln = time.time()
if bn256.g1_compress(M) == decrypted_message:
    print('ElGamal: ', telgln-telgl)
