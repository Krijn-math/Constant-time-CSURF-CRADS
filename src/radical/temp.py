# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 14:47:55 2021

@author: Krijn
"""

###TEMP
from random import randrange

e_2 = 150
fpadd = 0
fpsqr = 0
fpmul = 0

def set_zero_ops():

    global fpadd 
    global fpsqr
    global fpmul
    fpadd -= fpadd
    fpsqr -= fpsqr
    fpmul -= fpmul

def show_ops(label, a, b, flag):

    global fpadd 
    global fpsqr
    global fpmul

    print("| %s: %7dM + %7dS + %7da" % (label, fpmul, fpsqr, fpadd), end="\t")

    return None

def get_ops():
    
    global fpadd 
    global fpsqr
    global fpmul

    return [fpmul, fpsqr, fpadd]

# Modular addition
def fp_add(a, b):

    global fpadd
    fpadd += 1
    return (a + b) % p

# Modular substraction
def fp_sub(a, b):
    
    global fpadd
    fpadd += 1
    return (a - b) % p

# Modular multiplication
def fp_mul(a, b):
    
    global fpmul
    fpmul += 1
    return (a * b) % p

# Modular squaring
def fp_sqr(a):
    
    global fpsqr
    fpsqr += 1
    #print(fpsqr)
    return (a ** 2) % p

# constant-time swap
def fp_cswap(x, y, b):

    z = list([x, y])
    z = list(z[::(1 - 2*b)])
    return z[0], z[1]

# Modular exponentiation
def fp_exp(a, e):

    return pow(a, e, p)

# Modular inverse
def fp_inv(a):

	'''
	g, x, y = xgcd(a, p)
	#if g != 1:
	#	raise ValueError
	return x % p
	'''
	return fp_exp(a, p - 2)

def sign(x):
    if x < 0: return -1
    else: return 1
p = 1441439

inv_2 = fp_inv(2)       # 1 / 2

sq_exp = (p + 1)//4             #equiv to taking square root
quart_exp = (p + 1)//8          #if p suitable, equiv to taking fourth root
tri_exp = (2*p - 1)//3          #equiv to taking third root
quint_exp = (3*p-2)//5          #equiv to taking fifth root
sept_exp = (4*p-3)//7           #equiv to taking seventh root
novem_exp = (5*p - 4)//9        #equiv to taking ninth root
#KEPT THESE FOR POSSIBLE CHANGES

def Montgomery_to_Tate_nine_pro(A, xP, zP): 
    '''
    This is a projective version of the original Montgomery to Tate function
    specifically adapted to provide the right projective coordinates for the
    degree-nine projective radical isogeny. This assumes the point is given
    in projective coordinates.
'''


    rX = xP
    rZ = zP
    rZ2 = fp_sqr(rZ)
    rZ3 = fp_mul(rZ2, rZ)
    
    tmp = fp_mul(rZ,A)                  #rZ*A
    tmp = fp_add(tmp, rX)               #rZ*A + rX
    
    t0 = fp_mul(rX, tmp)                #rX*(rZ*A + rX)
    t0 = fp_add(t0, rZ2)
    t0 = fp_mul(rX, t0)
    t0 = fp_mul(rZ3, t0)                #rZ3*rX*(rZ2 + rX(rZA + rX))
    t0 = fp_exp(t0, sq_exp)
     
    twot = fp_add(t0,t0)                 #2*t0
    twot2 = fp_sqr(twot)                 #(2*t0)^2
    twot4 = fp_sqr(twot2)                #(2*t0)^4
    twot8 = fp_sqr(twot4)                #(2*t0)^8

    s0 = fp_add(tmp, tmp)               #2*rZ*A + 2*rX
    s0 = fp_add(s0, rX)
    s0 = fp_mul(rX, s0)
    s0 = fp_add(rZ2, s0)                #s0 = rZ^2 + 2*A*rX*rZ + 3*rX^2

    u0 = fp_add(rX, rX)
    u0 = fp_add(tmp, u0)
    u0 = fp_mul(twot2, u0)
    tmp2 = fp_sqr(s0)
    tmp2 = fp_mul(rZ3, tmp2)
    u0 = fp_sub(u0, tmp2)               #rZ*(2t0)^2*A + 3*(2t0)^2*rX - rZ^3*S0^2

    u3 = fp_sqr(u0)
    u3 = fp_mul(u3, u0)
    
    Mt = fp_mul(s0, u0)
    Mt = fp_add(Mt, Mt)                 #2*s0*u0
    
    tmp3 = fp_mul(rZ3, Mt)              #rZ3*Mt
    M = fp_sub(tmp3, twot4)
    M = fp_sqr(M)                       #(rZ^4*Mt - (2t0)^4)^2
    
    N = fp_mul(tmp3, twot4)
    N1 = fp_mul(rZ3, u3)
    N = fp_sub(N, N1)
    N = fp_sub(N, twot8)                #rZ^3 * (2t0)^4 * Mt - rZ^3*u0^3 - (2t0)^8
    
    return M, N

def Montgomery_to_Tate_nine_pro2(X, Z, xP, zP): 
    '''
    This is a projective version of the original Montgomery to Tate function
    specifically adapted to provide the right projective coordinates for the
    degree-nine projective radical isogeny. This assumes the point is given
    in projective coordinates.
'''

    rX = xP
    rZ = zP
    rX2 = fp_sqr(rX)
    rZ2 = fp_sqr(rZ)
    rZ3 = fp_mul(rZ2, rZ)
    
    Z2 = fp_sqr(Z)
    Z3 = fp_mul(Z2, Z)
    
    rXZ = fp_mul(rX, Z)
    
    tmp = fp_mul(rZ, X)
    tmp = fp_add(tmp, rXZ)              #rZ*A + rX
    
    tp = fp_mul(rXZ, tmp)                #rX*(rZ*A + rX)
    tp = fp_add(tp, fp_mul(Z2,rZ2))
    tp = fp_mul(rX, tp)
    tp = fp_mul(rZ3, tp)                #rZ3*rX*(rZ2 + rX(rZA + rX))
    tp = fp_mul(tp, Z2)
    tp = fp_exp(tp, sq_exp)
   
    twot = fp_add(tp,tp)                 #2*t0
    twot2 = fp_sqr(twot)                 #(2*t0)^2
    twot4 = fp_sqr(twot2)                #(2*t0)^4
    twot8 = fp_sqr(twot4)                #(2*t0)^8

    s1 = fp_mul(rZ2, Z)
    s2 = fp_mul(rX, rZ)
    s2 = fp_mul(X, s2)
    s2 = fp_add(s2, s2)
    s3 = fp_mul(rX2, Z)
    s3 = fp_add(s3, fp_add(s3, s3))
    sp = fp_add(s1, fp_add(s2, s3))

    u1 = fp_mul(rZ, X)
    u1 = fp_mul(u1, twot2)
    u2 = fp_mul(twot2, rXZ)
    u2 = fp_add(u2, fp_add(u2, u2))
    u3 = fp_sqr(sp)
    u3 = fp_mul(u3, rZ3)
    u3 = fp_mul(u3, Z3)
    up = fp_add(u1, u2)
    up = fp_sub(up, u3)

    u3 = fp_sqr(up)
    u3 = fp_mul(u3, up)
    
    Mtp = fp_mul(sp, up)
    Mtp = fp_add(Mtp, Mtp)                 #2*s0*u0

    Mp = fp_mul(rZ3, Mtp)              #rZ3*Mt
    Mp = fp_mul(Mp, Z2)
    Mp = fp_sub(Mp, twot4)
    Mp = fp_sqr(Mp)                       #(rZ^4*Mt - (2t0)^4)^2

    N1 = fp_mul(rZ3, Mtp)
    N1 = fp_mul(N1, twot4)
    N1 = fp_mul(N1, Z2)
    N2 = fp_mul(rZ3, u3)
    N2 = fp_mul(N2, Z)
    Np = fp_sub(N1, N2)
    Np = fp_sub(Np, twot8)
   
    return Mp, Np

B = randrange(p)
alf = randrange(p)
xP = randrange(p)
zP = randrange(p)

X, Z = Montgomery_to_Tate_nine_pro(B, xP, zP)
print(fp_mul(X, fp_inv(Z)))

X, Z = Montgomery_to_Tate_nine_pro2(fp_mul(alf, B), alf, xP, zP)
print(fp_mul(X, fp_inv(Z)))
