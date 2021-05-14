from framework import *
exec("from %s import *" % setting.algorithm)

print("// The exponents bounds and private keys go from the smallest to the largest small odd primes ell_i\'s\n")
if (setting.algorithm == 'csurf'):
    if setting.radicals:
        privkey_print('Exponents bounds', [e_2] + [e_3] + [e_5] + [e_7] + m[3:])
    else:
        privkey_print('Exponents bounds', [e_2] + m)
else:
    privkey_print('Exponents bounds', m)

print("\n// ===================== \033[0;33mPublic Key Generation\033[0m")

print("// --- \033[0;35mAlice\033[0m")
a, A, VELU_TIME, RAD_TIME, RUNNING_TIME = keygen()   # Alice
privkey_print('sk_a', a)
pubkey_print('pk_a', A)
print("// ~~~ Keygen Velu:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (VELU_TIME[0] / (10.0**6), VELU_TIME[1] / (10.0**6), VELU_TIME[2] / (10.0**6), measure(VELU_TIME) / (10.0**6)) )
print("// ~~~ Keygen Rads:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RAD_TIME[0] / (10.0**6), RAD_TIME[1] / (10.0**6), RAD_TIME[2] / (10.0**6), measure(RAD_TIME) / (10.0**6)) )
print("// ~~~ Keygen cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )


print("\n// --- \033[0;34mBob\033[0m")
b, B,  VELU_TIME, RAD_TIME, RUNNING_TIME = keygen()   # Bob
privkey_print('sk_b', b)
pubkey_print('pk_b', B)
print("// ~~~ Keygen Velu:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (VELU_TIME[0] / (10.0**6), VELU_TIME[1] / (10.0**6), VELU_TIME[2] / (10.0**6), measure(VELU_TIME) / (10.0**6)) )
print("// ~~~ Keygen Rads:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RAD_TIME[0] / (10.0**6), RAD_TIME[1] / (10.0**6), RAD_TIME[2] / (10.0**6), measure(RAD_TIME) / (10.0**6)) )
print("// ~~~ Keygen cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

print("\n// ===================== \033[0;33mSecret Sharing Computation\033[0m")
print("// --- \033[0;35mAlice\033[0m")
ss_a,  VELU_TIME, RAD_TIME, RUNNING_TIME = derive(a, B) # Alice
pubkey_print('ss_a', ss_a)
print("// ~~~ Derive Velu:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (VELU_TIME[0] / (10.0**6), VELU_TIME[1] / (10.0**6), VELU_TIME[2] / (10.0**6), measure(VELU_TIME) / (10.0**6)) )
print("// ~~~ Derive Rads:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RAD_TIME[0] / (10.0**6), RAD_TIME[1] / (10.0**6), RAD_TIME[2] / (10.0**6), measure(RAD_TIME) / (10.0**6)) )
print("// ~~~ Derive cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

print("\n// --- \033[0;34mBob\033[0m")
ss_b,  VELU_TIME, RAD_TIME, RUNNING_TIME = derive(b, A) # Bob
pubkey_print('ss_b', ss_b)
print("// ~~~ Derive Velu:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (VELU_TIME[0] / (10.0**6), VELU_TIME[1] / (10.0**6), VELU_TIME[2] / (10.0**6), measure(VELU_TIME) / (10.0**6)) )
print("// ~~~ Derive Rads:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RAD_TIME[0] / (10.0**6), RAD_TIME[1] / (10.0**6), RAD_TIME[2] / (10.0**6), measure(RAD_TIME) / (10.0**6)) )
print("// ~~~ Derive cost:\t%2.3fM + %2.3fS + %2.3fa = %2.3fM;" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

try:
    assert(ss_a == ss_b)
    print('\n\x1b[0;30;43m' + 'Successfully passed!' + '\x1b[0m')
except:
    raise TypeError('\x1b[0;30;41m' + 'Great Scott!... The sky is falling. NOT PASSED!!!' + '\x1b[0m')