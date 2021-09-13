from src.radical.crad_general import *
'''material needed for radical 4-isogenies'''
'''verified against Magma implementation of CRAD paper'''

#e_2 = 8 #274 #256

def act_with_two_on_Montgomery(A, exp):
    '''
    Given a Mont^+ coefficient, applies an isogeny to the surface, then applies 
    horizontal 2^exp-isogenies in the direction of the sign of exp. 
    Then applies an isogeny back to the floor. Due to the existence of 
    four-isogenies, this function is rarely used.
    '''
    
    A = Montgomery_to_Montgomery_min(A)
    A = (A * sign(exp)) % p         #
    
    A2 = fp_sqr(A)
    delta = fp_add(A2, 4)
    delta = crad_sq(delta)
    
    A_inv = fp_add(A, delta)
    A_inv = crad_inv(A_inv)           #1/(A + delta)
    
    A = fp_sub(A, delta)
    A = fp_sub(A, delta)
    A = fp_sub(A, delta)            #(A - 3*delta)
    A = fp_add(A,A)
    A = fp_mul(A, A_inv)            #2*(A - 3*delta)/(A + delta)
    
    for i in range(0, abs(exp) - 1):    #need to check VELU here
        eps = fp_sqr(A)
        eps = fp_sub(eps, 4)
        eps = crad_sq(eps)              #(A^2 - 4)^sq_exp
        Atmp = fp_sub(A, eps)
        Atmp = fp_mul(A, Atmp)
        Atmp = fp_sub(3, Atmp)
        A = fp_add(Atmp, Atmp)          #2*(3 - A*(A - eps))
    
    # Dummy isogenies (as used in MCR and OAYT-style)
    B = A
    for i in range( abs(exp)-1, e_2):
        eps = fp_sqr(B)
        eps = fp_sub(eps, 4)
        eps = crad_sq(eps)              #(A^2 - 4)^sq_exp
        Atmp = fp_sub(B, eps)
        Atmp = fp_mul(B, Atmp)
        Atmp = fp_sub(3, Atmp)
        B = fp_add(Atmp, Atmp)  

    eps = fp_sqr(A)
    eps = fp_sub(eps, 4)
    eps = crad_sq(eps)                 #(A^2 - 4)^sq_exp
    scalefac = fp_add(eps, A)
    scalefac = fp_mul(eps, scalefac)
    scalefac = fp_mul(scalefac, inv_2)
    scalefac = crad_sq(scalefac)       #(eps*(eps + A)/2)^sq_exp
    
    A_inv = fp_add(scalefac, scalefac)
    A_inv = crad_inv(A_inv)
    A = -A % p
    A = fp_sub(A, eps)
    A = fp_sub(A, eps)
    A = fp_sub(A, eps)
    A = fp_mul(A, A_inv)                #(-A - 3*eps)/(2*scalefac)
    
    A = (A * sign(exp)) % p         
    
    return Montgomery_min_to_Montgomery(A)

def alt_single(A):
    '''
    This function takes a Montgomery+ coefficient and applies a single
    2-isogeny.
    '''
    
    d = fp_add(A, 2)
    d = crad_sq(d)
    
    D = fp_sub(d, 2)
    D = fp_sqr(D)
    D = fp_sqr(D)
    
    C = fp_add(D, D)
    C = fp_add(C, C)
    
    inv = fp_sub(2, A)
    inv = fp_sqr(inv)
    inv = crad_inv(inv)
    
    out = fp_mul(C, inv)
    return fp_sub(2, out)

def alt_single_pro(X, Z):
    '''
    This function takes a Montgomery+ coefficient and applies a single
    2-isogeny, using projective coordinates
    '''
    
    Z2 = fp_sqr(Z)
    Z3 = fp_mul(Z, Z2)
    Z6 = fp_sqr(Z3)
    
    alfa = fp_add(Z, Z)
    alfa = fp_add(alfa, X)
    alfa = fp_mul(alfa, Z3)
    alfa = crad_sq(alfa) #((2Z + X)*Z^3)^sq_exp
        
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
    '''

    C = crad_quart(-A)                  #(-A)^quart_exp
    C4 = fp_add(C,C)                    #2*A^quart_exp
    C4 = fp_add(C4, C4)                 #4*A^quart_exp
    C42 = fp_mul(C4, C)                 #4*A^quart_exp*A^quart_exp
    output = fp_add(C42, 1)             #C42 + 1
    output = fp_mul(C, output)         #C * (C42 + 1)
    out_inv = fp_add(C42, C4)
    out_inv = fp_add(out_inv, 1)        #C42 - C4 + 1
    out_inv = fp_sqr(out_inv)           #(C42 - C4 + 1)^2
    out_inv = crad_inv(out_inv)           #1/((C42 - C4 + 1)^2)
    output = fp_mul(output, out_inv)    #-C*(C42+1)/((C42-C4+1)^2)
    
    return output

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

	alpha = crad_quart(-X)          # sqrt[4](-X)
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
    '''
	Degree-4 radical isogeny using projective coordinates.
	This is an adaption of the previous function, using the projective coordinates
    A  = (X:Z) and using the fact that (X:Z) = (XZ^7 : Z^8), thus taking the 4-throot
    alpha of -A in projective coordinates becomes (4-throot(XZ^7) : Z^2), where 
    the 4-th root is defined to be such that it takes the value a that is itself a
    quadratic residue
    '''
    Z2 = fp_sqr(Z)
    Z4 = fp_sqr(Z2)
    Z7 = fp_mul(Z4, fp_mul(Z2, Z))
    
    alpha = fp_mul(X, Z7)
    alpha = fp_sub(0, alpha)
    alpha = crad_quart(alpha)
           
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

def single_two_on_Montgomery_min(A):
    '''
    Applies a single isogeny of degree 2 on the surface
    '''
        
    A2 = fp_sqr(A)
    delta = fp_add(A2, 4)
    delta = crad_sq(delta)
    
    A_inv = fp_add(A, delta)
    A_inv = crad_inv(A_inv)           #1/(A + delta)
    
    A = fp_sub(A, delta)
    A = fp_sub(A, delta)
    A = fp_sub(A, delta)            #(A - 3*delta)
    A = fp_add(A,A)
    A = fp_mul(A, A_inv)            #2*(A - 3*delta)/(A + delta)
    
    eps = fp_sqr(A)
    eps = fp_sub(eps, 4)
    eps = crad_sq(eps)              #(A^2 - 4)^sq_exp
    Atmp = fp_sub(A, eps)
    Atmp = fp_mul(A, Atmp)
    Atmp = fp_sub(3, Atmp)
    A = fp_add(Atmp, Atmp)          #2*(3 - A*(A - eps))
    
    eps = fp_sqr(A)
    eps = fp_sub(eps, 4)
    eps = crad_sq(eps)                 #(A^2 - 4)^sq_exp
    scalefac = fp_add(eps, A)
    scalefac = fp_mul(eps, scalefac)
    scalefac = fp_mul(scalefac, inv_2)
    scalefac = crad_sq(scalefac)       #(eps*(eps + A)/2)^sq_exp
    
    A_inv = fp_add(scalefac, scalefac)
    A_inv = crad_inv(A_inv)
    A = -A % p
    A = fp_sub(A, eps)
    A = fp_sub(A, eps)
    A = fp_sub(A, eps)
    A = fp_mul(A, A_inv)                #(-A - 3*eps)/(2*scalefac)
        
    return (A)

def act_with_four_on_Montgomery_old(A, exp):
    '''
    Given a Mont^+ coefficient, applies a isogeny to the surface, then applies 
    horizontal 4^exp-isogenies = 2^(2*exp) in the direction of the sign of exp. 
    Then applies an isogeny back to the floor.
    
'''
    
    A = Montgomery_to_Montgomery_min(A) 

    A = (A * sign(exp)) % p         #
    d = fp_sqr(A)                   #A^2
    d = fp_add(d, 4)                #A^2 + 4
    d = crad_sq(d)                  #(A^2 + 4)^sq_exp
    D = fp_sub(A,d)                 #A-d
    D = fp_sqr(D)                   #(A-d)^2
    D = fp_add(D,4)                 #(A-d)^2 + 4
    D = crad_sq(D)                  #((A-d)^2 + 4)^sq_exp
    r = fp_add(D,d)                 #D + d
    r = fp_sub(r,A)                 #D + d - A
    twoinv = inv_2              #1/2
    r = fp_mul(r,twoinv)            #(D + d - A)/2
    t = fp_add(r, A)                #r + A
    t = fp_mul(t, r)                #(r + A)*r
    t = fp_sub(t, 1)                #(r + A)*r - 1
    t = fp_mul(t,r)                 #((r + A)*r - 1)*r
    t = crad_sq(t)                  #(((r + A)*r - 1)*r)^sq_exp

    _, A = Montgomery_min_to_Tate_four(A, r, t) #function returns M and N, need only N, M = 1 always holds
    
    X, Z = A, 1
    for i in range(0, abs(exp) // 2, 1):
        X, Z = four_iso_projective(X, Z)
        #X = four_iso(X)

    # Dummy isogenies (an easy patch, which must be modified in order to avoid dummy-isogenies)
    B1, B2 = X, Z
    for i in range(abs(exp) // 2, e_2 // 2, 1):        # Exponent bound of 2 is assumed to be e_2, thus c_2 // 2 degree-4 radical isogenies
        B1, B2  = four_iso_projective(B1, B2)
        #B = four_iso(X)

    # --- Recall, A = X / (Z^8)
    Z = crad_inv(Z)
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
    dum = single_two_on_Montgomery_min(output)
    output, dum = fp_cswap(output, dum, exp % 2)

    output = (output * sign(exp)) % p

    return Montgomery_min_to_Montgomery(output)

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
    d = crad_sq(d)                  #(A^2 + 4)^sq_exp    
    D = fp_sub(A,d)                 #A-d
    D = fp_sqr(D)                   #(A-d)^2
    D = fp_add(D,4)                 #(A-d)^2 + 4
    D = crad_sq(D)                  #((A-d)^2 + 4)^sq_exp
    r = fp_add(D,d)                 #D + d
    r = fp_sub(r,A)                 #D + d - A
    r = fp_mul(r,inv_2)            #(D + d - A)/2
    t = fp_add(r, A)                #r + A
    t = fp_mul(t, r)                #(r + A)*r
    t = fp_sub(t, 1)                #(r + A)*r - 1
    t = fp_mul(t,r)                 #((r + A)*r - 1)*r
    t = crad_sq(t)                  #(((r + A)*r - 1)*r)^sq_exp
    
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
    Z = crad_inv(Z)
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
    dum = alt_single(output)
    output, dum = fp_cswap(output, dum, exp % 2)

    output = (output * sign(exp)) % p

    return Montgomery_min_to_Montgomery(output)

def act_with_four_on_Montgomery_pro(AX, AZ, exp):
    '''
    Given a Mont^+ coefficient, applies a isogeny to the surface, then applies 
    horizontal 4^exp-isogenies = 2^(2*exp) in the direction of the sign of exp. 
    Then applies an isogeny back to the floor.
    
    Uses a specific version of moving from Mont- coefficients to the Tate 
    normal form, to save costs, as we only need a single output value of that
    function.
    
'''
    
    AX, AZ = Montgomery_to_Montgomery_min_pro(AX, AZ) 

    AX = (AX * sign(exp)) % p         #
        
    AZ2 = fp_sqr(AZ)
    _2AZ2 = fp_add(AZ2, AZ2)
    _4AZ2 = fp_add(_2AZ2, _2AZ2)
    
    dp = fp_sqr(AX)
    dp = fp_add(dp, _4AZ2)
    dp = crad_sq(dp)
     
    Dp = fp_sub(AX, dp)
    Dp = fp_sqr(Dp)
    Dp = fp_add(Dp, _4AZ2)
    Dp = crad_sq(Dp)
    
    rp = fp_add(Dp, dp)
    rp = fp_sub(rp, AX)
    rp = fp_mul(rp, inv_2)
    
    tp = fp_add(rp, AX)
    tp = fp_mul(tp, rp)
    tp = fp_sub(tp, AZ2)
    tp = fp_mul(tp, rp)
    tp = fp_mul(tp, AZ)
    tp = crad_sq(tp)    

    X, Z = Montgomery_min_to_Tate_four_pro(AX, AZ, rp, tp) #function only N, M = 1 always holds
        
    for i in range(0, abs(exp) // 2, 1):
        X, Z = four_iso_projective(X, Z)
        #X = four_iso(X)

    # Dummy isogenies (an easy patch, which must be modified in order to avoid dummy-isogenies)
    B1, B2 = X, Z
    for i in range(abs(exp) // 2, e_2 // 2, 1):        # Exponent bound of 2 is assumed to be e_2, thus c_2 // 2 degree-4 radical isogenies
        B1, B2  = four_iso_projective(B1, B2)
        #B = four_iso(X)
    
    X, Z = Tate_four_to_Montgomery_min_pro(X, Z)
    
    
    """
    if exp % 2:
        dum = single_two_on_Montgomery_min(A)
    else:
        A = single_two_on_Montgomery_min(A)
    """
    X, Z = Montgomery_min_to_Montgomery_pro(X, Z)
    
    dumX, dumZ = alt_single_pro(X, Z)
    X, dumX = fp_cswap(X, dumX, exp % 2)
    Z, dumZ = fp_cswap(Z, dumZ, exp % 2)


    X = (X * sign(exp)) % p
    
    out, inv = X, Z
    return fp_mul(out, crad_inv(inv))