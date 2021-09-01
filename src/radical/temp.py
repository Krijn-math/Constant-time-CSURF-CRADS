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

def single_two_on_Montgomery_min(A):
    '''
    Applies a single isogeny of degree 2 on the surface
    '''
        
    A2 = fp_sqr(A)
    delta = fp_add(A2, 4)
    delta = fp_exp(delta, sq_exp)
    
    A_inv = fp_add(A, delta)
    A_inv = fp_inv(A_inv)           #1/(A + delta)
    
    A = fp_sub(A, delta)
    A = fp_sub(A, delta)
    A = fp_sub(A, delta)            #(A - 3*delta)
    A = fp_add(A,A)
    A = fp_mul(A, A_inv)            #2*(A - 3*delta)/(A + delta)
    
    eps = fp_sqr(A)
    eps = fp_sub(eps, 4)
    eps = fp_exp(eps, sq_exp)              #(A^2 - 4)^sq_exp
    Atmp = fp_sub(A, eps)
    Atmp = fp_mul(A, Atmp)
    Atmp = fp_sub(3, Atmp)
    A = fp_add(Atmp, Atmp)          #2*(3 - A*(A - eps))
    
    eps = fp_sqr(A)
    eps = fp_sub(eps, 4)
    eps = fp_exp(eps, sq_exp)                 #(A^2 - 4)^sq_exp
    scalefac = fp_add(eps, A)
    scalefac = fp_mul(eps, scalefac)
    scalefac = fp_mul(scalefac, inv_2)
    scalefac = fp_exp(scalefac, sq_exp)       #(eps*(eps + A)/2)^sq_exp
    
    A_inv = fp_add(scalefac, scalefac)
    A_inv = fp_inv(A_inv)
    A = -A % p
    A = fp_sub(A, eps)
    A = fp_sub(A, eps)
    A = fp_sub(A, eps)
    A = fp_mul(A, A_inv)                #(-A - 3*eps)/(2*scalefac)
        
    return (A)

def alt_single(A):
    
    d = fp_add(A, 2)
    d = fp_exp(d, sq_exp)
    
    D = fp_sub(d, 2)
    D = fp_sqr(D)
    D = fp_sqr(D)
    
    C = fp_add(D, D)
    C = fp_add(C, C)
    
    inv = fp_sub(2, A)
    inv = fp_sqr(inv)
    inv = fp_inv(inv)
    
    out = fp_mul(C, inv)
    return fp_sub(2, out)

def alt_single_pro(X, Z):
    
    Z2 = fp_sqr(Z)
    Z3 = fp_mul(Z, Z2)
    Z6 = fp_sqr(Z3)
    
    alfa = fp_add(Z, Z)
    alfa = fp_add(alfa, X)
    alfa = fp_mul(alfa, Z3)
    alfa = fp_exp(alfa, sq_exp) #((2Z + X)*Z^3)^sq_exp
        
    inv = fp_add(Z, Z)
    inv = fp_sub(inv, X)
    inv = fp_sqr(inv)
    inv = fp_mul(inv, Z6)      #(2Z - X)^2*Z^6
    
    out = fp_sub(alfa, Z2)
    out = fp_sub(out, Z2)
    out = fp_sqr(out)
    out = fp_sqr(out)
    out = fp_add(out, out)
    tmp = fp_add(out, out)
    
    out = fp_add(inv, inv)
    out = fp_sub(out, tmp)      #2*inv - 4*(alfa - 2Z^2)^4
    
    return out, inv

def four_iso(A):
    '''
    Adapted from the Radical isogenies code: 
    Applies a single radical isogeny of degree 4 by mapping A to output.
    Due to the limitations of constant-time implementations, we obtain no speed-up
    with 4-isogenies, so we use 2-isogenies. 
    '''
    C = fp_exp(-A, quart_exp)                  #(-A)^quart_exp
    C4 = fp_add(C,C)                    #2*A^quart_exp
    C4 = fp_add(C4, C4)                 #4*A^quart_exp
    C42 = fp_mul(C4, C)                 #4*A^quart_exp*A^quart_exp
    output = fp_add(C42, 1)             #C42 + 1
    output = fp_mul(C, output)         #C * (C42 + 1)
    out_inv = fp_add(C42, C4)
    out_inv = fp_add(out_inv, 1)        #C42 + C4 + 1
    out_inv = fp_sqr(out_inv)           #(C42 + C4 + 1)^2

    out_inv = fp_inv(out_inv)           #1/((C42 + C4 + 1)^2)
    output = fp_mul(output, out_inv)    #-C*(C42+1)/((C42-C4+1)^2)
       
    return output
###TEMP

###NEW
def four_iso_projective_old(X, Z):
	'''
	Degree-4 radical isogeny using projective coordinates.
	In particular, we are writing A as X*Z^77 / Z^8 for some field elements X and Z. 
    The first isogeny is computed with X = A and Z = 1. 
    this representations allows you to re-use Z^2 as the fourth root of Z^8 from the
    previous iteration, saving one exponentiation in total.
	'''

	Z = fp_sqr(Z)					# Notice (Z^8) ^ [(p+1)//8] = Z ^ 2
	Z_squared = fp_sqr(Z)           # Z ^ 2

	alpha = fp_exp(-X, quart_exp)          # sqrt[4](-X)
	X_new = fp_add(alpha, alpha)    # 2 * alpha

	Z_new = fp_add(X_new, Z)            # 2 * alpha + Z

	X_new = fp_sqr(X_new)               #  4 * alpha ^ 2
	X_new = fp_add(X_new, Z_squared)    # (4 * alpha ^ 2 + Z ^ 2)
	X_new = fp_mul(X_new, alpha)        # (4 * alpha ^ 2 + Z ^ 2) * alpha
	X_new = fp_mul(X_new, Z)            # (4 * alpha ^ 2 + Z ^ 2) * alpha * Z

	tmp = fp_sqr(Z_new)                 # (2 * alpha + Z) ^ 2
	tmp = fp_sqr(tmp)                   # (2 * alpha + Z) ^ 4
	X_new = fp_mul(tmp, X_new)          # (2 * alpha + Z) ^ 4 * (4 * alpha ^ 2 + Z ^ 2) * alpha * Z

	# After all computations, we have The Tate normal form A-coefficient equals X_new / (Z_new ^ 8)
	# Thus, the ouput can be taken as input in an iterative chain of degree-4 radical isogenies
	return X_new, Z_new

def four_iso_projective(X, Z):
      
    Z2 = fp_sqr(Z)
    Z4 = fp_sqr(Z2)
    Z7 = fp_mul(Z4, fp_mul(Z2, Z))
    
    alpha = fp_mul(X, Z7)
    alpha = fp_sub(0, alpha)
    alpha = fp_exp(alpha, quart_exp)
           
    alpha2 = fp_sqr(alpha)
    out = fp_add(alpha2, alpha2)
    out = fp_add(out, out)
    tmp = fp_add(out, Z4)               #4*a^2 + Z^4
    out = fp_mul(tmp, alpha)  
    out = fp_mul(out, Z2)
    
    inv = fp_mul(alpha, Z2)
    inv = fp_add(inv, inv)
    inv = fp_add(inv, inv)
    inv = fp_add(inv, tmp)
    inv = fp_sqr(inv)                   #(4*a^2 + Z^4 + 4*a*z^2)^2

    return out, inv

def Montgomery_min_to_Montgomery_pro(AX, AZ):  
    ''' From the Radical isogenies code: 
    This function transforms
        FROM: a Montgomery^- coefficient A on the surface
                    y^2 = x^3 + A*x^2 - x.
        TO: a Montgomery coefficient A on the floor
                    y^2 = x^3 + A*x^2 + x
    obtained from the 2-isogeny with kernel <(0,0)>.
'''
    output = fp_add(-AX,-AX)     
    AX2 = fp_sqr(AX)              
    AZ2 = fp_sqr(AZ)
    inv = fp_add(AZ2, AZ2)
    inv = fp_add(inv, inv)
    inv = fp_add(inv, AX2)         
    inv = fp_exp(inv, sq_exp)      

    return output, inv

def Montgomery_to_Montgomery_min_pro(AX, AZ):
    '''
    From the Radical isogenies code:   
    This function transforms
        FROM: a Montgomery coefficient A on the floor
                    y^2 = x^3 + A*x^2 + x
        TO: a Montgomery^- coefficient A on the surface
                    y^2 = x^3 + A*x^2 - x.
    obtained from the 2-isogeny with kernel <(0,0)>. 
'''
  
    output = fp_add(-AX,-AX)
    output = fp_mul(output, AZ)
    
    AX2 = fp_sqr(AX)               
    AZ2 = fp_sqr(AZ)
    inv = fp_add(AZ2, AZ2)
    inv = fp_add(inv, inv)
    inv = fp_sub(inv, AX2)
    inv = fp_mul(inv, AZ2)
    inv = fp_exp(inv, sq_exp)      
    
    return output, inv

def Montgomery_min_to_Tate_four_pro(AX, AZ, rp, tp): 
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
    
    our function only returns N, in projective coordinates, to save costs
    
    notice r = rp/Z, t = tp/Z^2 
'''
    AZ2 = fp_sqr(AZ)
    AZ4 = fp_sqr(AZ2)

    tau = fp_add(tp, tp)
    
    sp = rp
    sp = fp_add(sp, fp_add(sp, sp))
    tmp = fp_add(sp, AX)
    sp = fp_add(tmp, AX)
    sp = fp_mul(rp, sp)
    sp = fp_sub(sp, AZ2)                 #-Z^2 + rp(2X + 3rp)
     
    up = fp_mul(fp_sqr(tau), tmp)
    up = fp_sub(up, fp_mul(AZ, fp_sqr(sp)))         #tau^2(X+3rp) - sp^2*Z
      
    NX = fp_mul(up, fp_sqr(up))         
    NX = fp_mul(NX, AZ)                    #up^3*Z^4
    
    tau4 = fp_sqr(fp_sqr(tau))
    NZ = fp_sqr(tau4)                     #tau^8

    return NX, NZ


def Tate_four_to_Montgomery_min_pro(AX, AZ):
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
    
    _2AX = fp_add(AX, AX)
    _4AX = fp_add(_2AX, _2AX)
    _4AX2 = fp_sqr(_4AX)
    _4AX3 = fp_mul(_4AX, _4AX2)
    _4AX6 = fp_sqr(_4AX3)
    _8AX = fp_add(_4AX, _4AX)

    Bp = fp_sub(AZ, _8AX)                   #Z - 8X
    
    eps = fp_add(Bp, _8AX)
    tmp = fp_sub(Bp, _8AX)
    eps = fp_mul(eps, tmp)
    fix = fp_sqr(_4AX)
    eps = fp_mul(eps, fix)
    eps = fp_exp(eps, sq_exp)               #((Bp + 8*X)(Bp - 8*X)*(4X)^2)^sq_exp

    _3eps = fp_add(eps, fp_add(eps, eps))
    ttmp = fp_mul(Bp, _4AX)
    output = fp_sub(_3eps, ttmp)            #-Bp*(4x) + 3*eps

    inv = fp_mul(Bp, _4AX)
    inv = fp_sub(inv, eps)
    inv = fp_mul(eps, inv)
    inv = fp_add(inv, inv)
    inv = fp_exp(inv, sq_exp)               #(2*eps*(BP(4X) - eps))^sq_exp
    
    return output, inv

###TODO

def act_with_four_on_Montgomery_pro(AX, AZ, exp):
    '''
    Given a Mont^+ coefficient, applies a isogeny to the surface, then applies 
    horizontal 4^exp-isogenies = 2^(2*exp) in the direction of the sign of exp. 
    Then applies an isogeny back to the floor.
    
'''
    
    AX, AZ = Montgomery_to_Montgomery_min_pro(AX, AZ) 

    AX = (AX * sign(exp)) % p         #
        
    AZ2 = fp_sqr(AZ)
    _2AZ2 = fp_add(AZ2, AZ2)
    _4AZ2 = fp_add(_2AZ2, _2AZ2)
    
    dp = fp_sqr(AX)
    dp = fp_add(dp, _4AZ2)
    dp = fp_exp(dp, sq_exp)
     
    Dp = fp_sub(AX, dp)
    Dp = fp_sqr(Dp)
    Dp = fp_add(Dp, _4AZ2)
    Dp = fp_exp(Dp, sq_exp)
    
    rp = fp_add(Dp, dp)
    rp = fp_sub(rp, AX)
    rp = fp_mul(rp, inv_2)
    
    tp = fp_add(rp, AX)
    tp = fp_mul(tp, rp)
    tp = fp_sub(tp, AZ2)
    tp = fp_mul(tp, rp)
    tp = fp_mul(tp, AZ)
    tp = fp_exp(tp, sq_exp)    

    X, Z = Montgomery_min_to_Tate_four_pro(AX, AZ, rp, tp) #function only N, M = 1 always holds
        
    for i in range(0, abs(exp) // 2, 1):
        X, Z = four_iso_projective(X, Z)
        #X = four_iso(X)

    # Dummy isogenies (an easy patch, which must be modified in order to avoid dummy-isogenies)
    B1, B2 = X, Z
    for i in range(abs(exp) // 2, e_2 // 2, 1):        # Exponent bound of 2 is assumed to be e_2, thus c_2 // 2 degree-4 radical isogenies
        B1, B2  = four_iso_projective(B1, B2)
        #B = four_iso(X)

    A = fp_mul(X, fp_inv(Z))
    
    X, Z = Tate_four_to_Montgomery_min_pro(X, Z)
    
    
    """
    if exp % 2:
        dum = single_two_on_Montgomery_min(A)
    else:
        A = single_two_on_Montgomery_min(A)
    """
    #RESTORE THESE
    dumX, dumZ = alt_single_pro(X, Z)
    X, dumX = fp_cswap(X, dumX, exp % 2)
    Z, dumZ = fp_cswap(Z, dumZ, exp % 2)


    X = (X * sign(exp)) % p
    
    return Montgomery_min_to_Montgomery_pro(X, Z)

def act_with_four_on_Montgomery(A, exp):
    '''
    Given a Mont^+ coefficient, applies a isogeny to the surface, then applies 
    horizontal 4^exp-isogenies = 2^(2*exp) in the direction of the sign of exp. 
    Then applies an isogeny back to the floor.
    
'''
    
    A = Montgomery_to_Montgomery_min(A) 

    A = (A * sign(exp)) % p         #
    d = fp_sqr(A)                   #A^2
    d = fp_add(d, 4)                #A^2 + 4
    d = fp_exp(d, sq_exp)                  #(A^2 + 4)^sq_exp    
    D = fp_sub(A,d)                 #A-d
    D = fp_sqr(D)                   #(A-d)^2
    D = fp_add(D,4)                 #(A-d)^2 + 4
    D = fp_exp(D, sq_exp)                  #((A-d)^2 + 4)^sq_exp
    r = fp_add(D,d)                 #D + d
    r = fp_sub(r,A)                 #D + d - A
    r = fp_mul(r,inv_2)            #(D + d - A)/2
    t = fp_add(r, A)                #r + A
    t = fp_mul(t, r)                #(r + A)*r
    t = fp_sub(t, 1)                #(r + A)*r - 1
    t = fp_mul(t,r)                 #((r + A)*r - 1)*r
    t = fp_exp(t, sq_exp)                  #(((r + A)*r - 1)*r)^sq_exp
    
    _, A = Montgomery_min_to_Tate_four(A, r, t) #function returns M and N, need only N, M = 1 always holds
    X, Z = A, 1
    
    for i in range(0, abs(exp) // 2, 1):
        X, Z = four_iso_projective_old(X, Z)
        #X = four_iso(X)

    # Dummy isogenies (an easy patch, which must be modified in order to avoid dummy-isogenies)
    B1, B2 = X, Z
    for i in range(abs(exp) // 2, e_2 // 2, 1):        # Exponent bound of 2 is assumed to be e_2, thus c_2 // 2 degree-4 radical isogenies
        B1, B2  = four_iso_projective_old(B1, B2)
        #B = four_iso(X)

    # --- Recall, A = X / (Z^8)
    Z = fp_inv(Z)
    Z = fp_sqr(Z)
    Z = fp_sqr(Z)
    Z = fp_sqr(Z)
    A = fp_mul(X, Z)

    output = Tate_four_to_Montgomery_min(A)
    """
    if exp % 2:
        dum = single_two_on_Montgomery_min(A)
    else:
        A = single_two_on_Montgomery_min(A)
    """
    #RESTORE THESE
    dum = alt_single(output)
    output, dum = fp_cswap(output, dum, exp % 2)

    output = (output * sign(exp)) % p

    return Montgomery_min_to_Montgomery(output)

#A = randrange(p)
#print(act_with_four_on_Montgomery(A, 3))



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

