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

#KEPT THESE FOR POSSIBLE CHANGES

def Montgomery_min_to_Montgomery(A):  
    ''' From the Radical isogenies code: 
    This function transforms
        FROM: a Montgomery^- coefficient A on the surface
                    y^2 = x^3 + A*x^2 - x.
        TO: a Montgomery coefficient A on the floor
                    y^2 = x^3 + A*x^2 + x
    obtained from the 2-isogeny with kernel <(0,0)>.
'''

    output = fp_add(-A,-A)        #(-2*A)
    inv = fp_sqr(A)               #(A^2)
    inv = fp_add(inv, 4)          #(A^2 + 4)
    inv = fp_exp(inv, sq_exp)            #(A^2 + 4)^sq_exp
    inv = fp_inv(inv)             #1/(A^2 + 4)^sq_exp
    output = fp_mul(output,inv)   #(-2*A)/((A^2 + 4)^sq_exp)
    return output

def Montgomery_to_Montgomery_min(A):
    '''
    From the Radical isogenies code:   
    This function transforms
        FROM: a Montgomery coefficient A on the floor
                    y^2 = x^3 + A*x^2 + x
        TO: a Montgomery^- coefficient A on the surface
                    y^2 = x^3 + A*x^2 - x.
    obtained from the 2-isogeny with kernel <(0,0)>. 
'''
  
    output = fp_add(-A,-A)        #(-2*A)
    inv = fp_sqr(A)               #(A^2)
    inv = fp_sub(4, inv)          #(4 - A^2)
    inv = fp_exp(inv, sq_exp)     #(4 - A^2)^sq_exp
    inv = fp_inv(inv)             #1/((4 - A^2)^sq_exp)
    output = fp_mul(output,inv)   #(-2*A)/((4 - A^2)^sq_exp)
    return output

def Montgomery_min_to_Tate_four(A, r, t): 
    '''
    From Radical Isogenies:
    This functions transforms
        FROM: Montgomery^- coefficient A
                y^2 = x^3 + A*x^2 - x
        TO: the coefficients of a Tate normal form of X1(4)
                y^2 + M*xy + N*y = x^3 + N*x^2
    
    In this case for X1(4), M = 1 will always hold.
    Note that [r,t] is the Fp-rational 4-torsion point such that the isogeny with kernel
    2*[r,t] maps [r,t] to a 2-torsion point that has Fp-rational halves as well.
'''
    
    tmp = fp_add(A, r)  
    twot_inv = fp_add(t,t)              #2*t
    twot_inv = fp_inv(twot_inv)         #1/(2*t)
    s = fp_add(tmp, tmp)                #2*A + 2*r
    s = fp_add(s, r)                    #2*A + 3*r
    s = fp_mul(r,s)                     #r*(2*A + 3*r)
    s = fp_add(-1,s)                    #-1 + r*(2*A + 3*r) TODO: CHECK WHY -1 HERE
    s = fp_mul(s, twot_inv)             #(-1 + r*(2*A + 3*r))/(2*t)

    u_inv = fp_add(r,r)                 #2*r
    u_inv = fp_add(tmp, u_inv)          #A + 3*r
    s_sq = fp_sqr(s)                    #s^2
    u_inv = fp_sub(u_inv, s_sq)         #A + 3*r - s^2
    u_inv = fp_mul(u_inv, twot_inv)     #((A + 3*r - s^2)/(2*t))
        
    M = fp_mul(s, u_inv)                #s*u_inv
    M = fp_add(M,M)                     #2*s*u_inv
    N = fp_mul(u_inv, u_inv)
    N = fp_mul(u_inv, N)                #u_inv^3
    N = fp_mul(t, N)                    #t*u_inv^3
    N = fp_add(N,N)                     #2*t*u_inv^3

    return M, N

def Tate_four_to_Montgomery_min(A):
    '''
    From Radical Isogenies:
    This function transforms
        FROM: coefficient A representing a Tate normal form 
                    y^2 + x*y + A*y = x^3 + A*x^2 
        TO: Montgomery^- coefficient A
                y^2 = x^3 + A*x^2 - x
    B is a Montgomery coefficient on the surface where we simply translated 2*(0,0) = (-A,0) to (0,0).
    The rest is a classical rescaling to obtain a Montgomery^- coefficient (just as in CSURF).
'''
    
    B = fp_add(A,A)                         #2*A
    B = fp_add(B,B)                         #4*A
    B = fp_inv(B)                           #1/4*A
    B = fp_sub(B,2)                         #1/4*A - 2
    
    eps = fp_sqr(B)                         #B^2
    eps = fp_sub(eps, 4)                    #B^2 - 4
    eps = fp_exp(eps, sq_exp)               #(B^2 - 4)^sq_exp
    
    output = fp_mul(3,eps)                  #3*eps
    output = fp_add(-B, output)             #(-B + 3*eps)
    out_inv = fp_sub(B,eps)                 #(B-eps)
    out_inv = fp_mul(eps, out_inv)          #eps*(B-eps)
    out_inv = fp_add(out_inv,out_inv)       #eps*(B-eps)*2
    out_inv = fp_exp(out_inv, sq_exp)       #(eps*(B-eps)*2)^sq_exp

    out_inv = fp_inv(out_inv)             #1/(#(eps*(B-eps)*2)^sq_exp)
    output = fp_mul(output, out_inv)        #(-B + 3*eps)/(#(eps*(B-eps)*2)^sq_exp)
    
    return output

