from src.radical.crad_general import *
from src.montgomery import isinfinity

'''material needed for radical 5-isogenies'''
#e_5 = 0

def five_iso_projective(X,Z):
    '''
	Degree-5 radical isogeny using projective coordinates.
	In particular, we are writing A as X*Z^4 / Z^5 for some field elements X and Z. 
    The first isogeny is computed with X = A and Z = 1. 
    this representations allows you to re-use Z as the fifth root of Z^5 from the
    previous iteration, saving one exponentiation in total.
    '''
    
    X = crad_quint(X)
    Z = Z
    
    X2 = fp_sqr(X)
    X3 = fp_mul(X, X2)
    X4 = fp_sqr(X2)
    
    Z2 = fp_sqr(Z)
    Z3 = fp_mul(Z, Z2)
    Z4 = fp_sqr(Z2)
    
    X3Z = fp_mul(X3, Z)
    X2Z2 = fp_mul(X2, Z2)
    XZ3 = fp_mul(X, Z3)
    
    num = fp_add(X4, fp_mul( 3, X3Z))
    num = fp_add(num, fp_mul( 4, X2Z2))
    num = fp_add(num, fp_mul( 2, XZ3))
    num = fp_add(num, Z4)
    num = fp_mul(X, num)                #X(X^4 + 3*X^3*Z + 4*X^2*Z^2 + 2*X*Z^3 + Z^4)
        
    denom = fp_add(X4, fp_mul( -2, X3Z))
    denom = fp_add(denom, fp_mul( 4, X2Z2))
    denom = fp_add(denom, fp_mul( -3, XZ3))
    denom = fp_add(denom, Z4)
    denom = fp_mul(Z, denom)            #Z(X^4 - 2*X^3*Z + 4*X^2*Z^2 - 3*X*Z^3 + Z^4)    
    
    denom2 = fp_sqr(denom)
    denom4 = fp_sqr(denom2)
    
    num = fp_mul(num, denom4)

    return num, denom
    
    
def pro_to_aff_five(X, Z):
    '''
    turns the projective output from five_iso_projective back to affine
    coordinates by computing A = X/Z^5
    '''
    Z2 = fp_sqr(Z)
    Z4 = fp_sqr(Z2)
    Z5 = fp_mul(Z, Z4)

    A = crad_inv(Z5)
    A = fp_mul(X,A)
    return A


def act_with_five_on_Montgomery(A, exp, Tp_proj = None):
    '''
    Adapted from the Radical isogenies code.
    Applies exp number of 5 isogenies to the Montgomery curve E_A on the floor.
    The sign of exp indicates the direction of the isogenies.
    We use the projective version of the 5-isogeny to save inversions.
    '''
    A = (A * sign(exp)) % p         #
    if Tp_proj == None:
        Tp_proj, _ = sampling_ell_order_point(A, 1)
        while isinfinity(Tp_proj):
            Tp_proj, _ = sampling_ell_order_point(A, 1)
        
        sign_exp = 1
    else:
        sign_exp = sign(exp)

    # Affine x-coordinate Tp
    xTp = crad_inv(Tp_proj[1])
    xTp = fp_mul(xTp, Tp_proj[0])
    xTp = (xTp * sign_exp) % p         # <--- isomorphism from A to -A with zera-trace points on A
    
    _, A = Montgomery_to_Tate(A, xTp)
    A = -A % p                      #needed for input five_iso
    
    X, Z = A, 1
    for i in range(0, abs(exp), 1):
        X, Z = five_iso_projective(X, Z)

    for i in range(abs(exp), e_5, 1):        # Exponent bound of 5 is assumed to be e_5, thus c_5 degree-5 radical isogenies
        B1, B2 = five_iso_projective(X, Z)
    
    A = pro_to_aff_five(X, Z)
    
    A = -A % p
    A2 = fp_sqr(A)
    A3 = fp_mul(A, A2)
    A4 = fp_sqr(A2)
    A5 = fp_mul(A, A4)
    
    a1 = fp_add(A,1)

    a2 = A
    a3 = a2
    
    #we do not apply a last Velu in these last coordinates
    a4 = 0
    a6 = 0
    
    # a4 = fp_mul(5, A3)                  #5*A^3
    # a4 = fp_add(a4, fp_mul(-10, A2))    #5*A^3 - 10*A^2
    # a4 = fp_add(a4, fp_mul( -5,  A))    #5*A^3 - 10*A^2 - 5*A
    
    # a6 = fp_add(A5, fp_mul(-10, A4))    #A^5 - 10A^4
    # a6 = fp_add(a6, fp_mul( -5, A3))    #A^5 - 10A^4 - 5A^3
    # a6 = fp_add(a6, fp_mul(-15, A2))    #A^5 - 10A^4 - 5A^3 - 15A^2
    # a6 = fp_sub(a6, A)                  #A^5 - 10A^4 - 5A^3 - 15A^2 - A
  
    output = Weier_to_Montgomery([a1, a2, a3, a4, a6])
    output = (output * sign(exp)) % p         #

    return output