from src.radical.crad_general import *
from src.montgomery import isinfinity

'''material needed for radical 7-isogenies'''

#e_7 = 0

def seven_iso(A):
    '''
    Adapted from the Radical isogenies code: 
    Applies a single radical isogeny of degree 7
    The input A represents an element of X1(7) in Tate normal form with
    (0,0) as 7-torsion point, 
    The output represents the Tate normal form that is obtained from the isogeny with kernel <(0,0)>.
    Due to the inversion required, the projective version of this isogeny
    is faster and is therefore used.
    '''
    
    A2 = fp_sqr(A)
    A3 = fp_mul(A2,A)
    A4 = fp_sqr(A2)
    A5 = fp_mul(A4, A)
    
    inv = fp_add(A5, fp_mul( -8, A4))       #A^5 - 8*A^4
    inv = fp_add(inv, fp_mul( 5, A3))       #A^5 - 8*A^4 + 5*A^3
    inv = fp_add(inv, A2)                   #A^5 - 8*A^4 + 5*A^3 + A^2
    inv = crad_inv(inv)                     #1/(#A^5 - 8*A^4 + 5*A^3 + A^2)
    
    C = fp_sub(A5, A4)
    C = crad_sept(C)                        #sqrt[7](A^5 - A^4)

    S1 = fp_add(C, 1)
    S1 = fp_mul(-7, S1)                     #-7C -7
    
    S2 = fp_mul(4, A2)
    S2 = fp_add(S2, fp_mul(-12, A))
    S2 = fp_add(S2, 2)                      #4*A^2 - 12*A + 2
    
    S3 = fp_add(A3, fp_mul(-10, A2))
    S3 = fp_add(S3, fp_mul(  4,  A))        #A^3 - 10A^2 + 4A
        
    S4 = fp_mul(-5, A3)
    S4 = fp_add(S4, A2)
    S4 = fp_add(S4, A)                      #-5A^3 + A2 + A
    
    S5 = fp_mul(3, A4)
    S5 = fp_add(S5, fp_mul(-9, A3))
    S5 = fp_add(S5, fp_mul( 5, A2))         #3A^4 - 9A^3 + 5A^2
    
    S6 = fp_mul(3, A5)
    S6 = fp_add(S6, fp_mul(-16, A4))
    S6 = fp_add(S6, fp_mul( 12, A3))         #3A^5 - 16A^4 + 12A^3
    
    output = fp_mul(S1, C)
    output = fp_add(output, S2)
    output = fp_mul(output, C)
    output = fp_add(output, S3)
    output = fp_mul(output, C)
    output = fp_add(output, S4)
    output = fp_mul(output, C)
    output = fp_add(output, S5)
    output = fp_mul(output, C)
    output = fp_add(output, S6)
    output = fp_mul(output, inv)

    return output

def seven_iso_projective(X, Z):
    '''
	Degree-7 radical isogeny using projective coordinates.
	In particular, we are writing A as X*Z^6 / Z^7 for some field elements X and Z. 
    The first isogeny is computed with X = A and Z = 1. 
    this representations allows you to re-use Z as the seventh root of Z^7 from the
    previous iteration, saving one exponentiation in total. Here, we need to
    calculate beta = sqrt[7](X^4*Z^2*(X-Z)).    
    the rest of the formulas are projective versions of the above
    '''
    
    X2 = fp_sqr(X)
    X3 = fp_mul(X2, X)
    X4 = fp_sqr(X2)
    X5 = fp_mul(X4, X)
    
    Z2 = fp_sqr(Z)
    Z3 = fp_mul(Z2, Z)
    
    X4Z = fp_mul(X4, Z)
    X3Z2 = fp_mul(X3, Z2)
    X2Z3 = fp_mul(X2, Z3)
    
    denom = fp_add(X5, fp_mul(-8, X4Z))
    denom = fp_add(denom, fp_mul(5, X3Z2))
    denom = fp_add(denom, X2Z3)
    denom = fp_mul(Z, denom)                #Z*(X^5 - 8*X^4*Z + 5*X^3*Z^2 + X^2*Z^3)
    
    XZ = fp_mul(X,Z)
    X2Z = fp_mul(X2, Z)
    X3Z = fp_mul(X3, Z)
    X4Z = fp_mul(X4, Z)
    
    XZ2 = fp_mul(X, Z2)
    X2Z2 = fp_mul(X2, Z2)
    
    beta = fp_sub(X, Z)
    beta = fp_mul(Z2, beta)
    beta = fp_mul(X4, beta)
    beta = crad_sept(beta)                   #sqrt[7](X^4*Z^2*(X-Z))
   
    S1 = fp_add(beta, Z)                     #(beta + Z)
    S1 = fp_mul(-7, S1)                      #-7*beta - 7*Z

    S2 = fp_mul(4, X2)
    S2 = fp_add(S2, fp_mul(-12, XZ))
    S2 = fp_add(S2, fp_mul(  2, Z2))         #(4X^2 - 12XZ + 2Z^2)

    S3 = fp_add(X3, fp_mul(-10, X2Z))
    S3 = fp_add(S3, fp_mul(  4, XZ2))       #(X^3 - 10X^2Z + 4XZ^2)
    
    S4 = fp_mul(-5, X3)
    S4 = fp_add(S4, X2Z)
    S4 = fp_add(S4, XZ2)
    S4 = fp_mul(S4, Z)                      #(-5X^3 + X^2Z + XZ^2)*Z
    
    S5 = fp_mul(3, X4)
    S5 = fp_add(S5, fp_mul(-9, X3Z))
    S5 = fp_add(S5, fp_mul( 5, X2Z2))
    S5 = fp_mul(S5, Z)                      #(3X^4 - 9X^3Z + 5X^2Z^2)*Z
    
    S6 = fp_mul(3, X5)
    S6 = fp_add(S6, fp_mul(-16, X4Z))
    S6 = fp_add(S6, fp_mul( 12, X3Z2))
    S6 = fp_mul(S6, Z)                      #(3X^5 - 16X^4Z + 12X^3Z^2)*Z
    
    num = fp_mul(S1, beta)
    num = fp_add(num, S2)
    num = fp_mul(num, beta)
    num = fp_add(num, S3)
    num = fp_mul(num, beta)
    num = fp_add(num, S4)
    num = fp_mul(num, beta)
    num = fp_add(num, S5)
    num = fp_mul(num, beta)
    num = fp_add(num, S6)
        
    return num, denom

def pro_to_aff_seven(X,Z):
    '''
    A  is represented as XZ^6/Z^7, so this turns it back to affine coordinates
    '''
    A = crad_inv(Z)
    A = fp_mul(X, A)
    
    return A

def act_with_seven_on_Montgomery(A, exp):
    '''
    Adapted from the Radical isogenies code.
    Applies exp number of 7 isogenies to the Montgomery curve E_A on the floor.
    The sign of exp indicates the direction of the isogenies.
    We use the projective version of the 7-isogeny to save inversions
    '''

    A = (A * sign(exp)) % p         #

    Tp_proj, _ = sampling_ell_order_point(A, 2)
    while isinfinity(Tp_proj):
        Tp_proj, _ = sampling_ell_order_point(A, 2)

    # Affine x-coordinate Tp
    xTp = crad_inv(Tp_proj[1])
    xTp = fp_mul(xTp, Tp_proj[0])
        
    M, N = Montgomery_to_Tate(A, xTp)
    #invA = fp_sub(M, 1)
    #invA = fp_inv(invA)
    #A = fp_mul(N, invA)                         #these steps are needed for affine isogenies
    
    
    X, Z = N, fp_sub(M,1)
    
    for i in range(0, abs(exp)):
        X, Z = seven_iso_projective(X, Z)

    for i in range(abs(exp), e_7, 1):        # Exponent bound of 7 is assumed to be e_7, thus c_7 degree-7 radical isogenies
        B1, B2  = seven_iso_projective(X, Z)
        
    A = pro_to_aff_seven(X, Z)
    
    A2 = fp_sqr(A)
    A3 = fp_mul(A, A2)
    A4 = fp_sqr(A2)
    A5 = fp_mul(A4, A)
    A6 = fp_sqr(A3)
    A7 = fp_mul(A6, A)
    A8 = fp_sqr(A4)
    A9 = fp_mul(A8, A)
    A10 = fp_sqr(A5)
    A11 = fp_mul(A10, A)
    
    a1 = fp_add(1, A)
    a1 = fp_sub(a1, A2)                 #1+A-A2
    
    a2 = fp_sub(A2, A3)                 #A2 - A3
    a3 = a2
    
    a4 = 0                              #we do not apply Velu in the last step
    
    a6 = 0
    
    # a4 = fp_mul(-5, A7)
    # a4 = fp_add(a4, fp_mul( 35, A5))
    # a4 = fp_add(a4, fp_mul(-70, A4))
    # a4 = fp_add(a4, fp_mul( 70, A3))
    # a4 = fp_add(a4, fp_mul(-35, A2))
    # a4 = fp_add(a4, fp_mul(  5,  A))    #(-5*A7 + 35*A5 - 70*A4 + 70*A3 - 35*A2 + 5*A)

    # a6 = fp_sub(A, A11)
    # a6 = fp_add(a6, fp_mul(  -8, A10))
    # a6 = fp_add(a6, fp_mul(  46,  A9))
    # a6 = fp_add(a6, fp_mul(-107,  A8))
    # a6 = fp_add(a6, fp_mul( 202,  A7))
    # a6 = fp_add(a6, fp_mul(-343,  A6))
    # a6 = fp_add(a6, fp_mul( 393,  A5))
    # a6 = fp_add(a6, fp_mul(-258,  A4))
    # a6 = fp_add(a6, fp_mul(  94,  A3))
    # a6 = fp_add(a6, fp_mul( -19,  A2))  #(-A11 - 8*A10 + 46*A9 - 107*A8 + 202*A7 - 343*A6 + 393*A5 - 258*A4 + 94*A3 - 19*A2 + A)

    output = Weier_to_Montgomery([a1, a2, a3, a4, a6])
    output = (output * sign(exp)) % p         #

    return output