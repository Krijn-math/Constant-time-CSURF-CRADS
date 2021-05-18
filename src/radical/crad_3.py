from src.radical.crad_general import *
from src.montgomery import isinfinity

'''material needed for radical 3 and 9-isogenies'''
'''computations verified against Magma implementation CRAD'''

#e_3 = 8

def Montgomery_to_Tate_three(A, xP):
    '''
    Adapted from the Radical isogenies code: 
    This functions transforms
        FROM: a Montgomery coefficient A on the floor
                    y^2 = x^3 + A*x^2 + x
        TO: the coefficients of a Tate normal form of X1(4)
                y^2 + M*xy + N*y = x^3 + N*x^2
    obtained by an isomorphism that puts a point with x-coordinate xP at (0,0), 
    where xP is assumed to be a 3-torsion point's x-coordinate. 
    The choice of sign for the y-coordinate corresponding to xP is arbitrary but fixed. 
    Using notation from the Tate normal form, M = 1-c and N = -b.
'''

    r = xP
    tmp = fp_add(A,r)                   #A+r
    t = fp_mul(r,tmp)                   #r*(A+r)
    t = fp_add(1, t)                    #(1 + r*(A+r))
    t = fp_mul(r,t)                     #r * (1 + r*(A+r))
    t = crad_sq(t)                      #(r*(1+r*(A+r)))^sq_exp
    twot_inv = fp_add(t,t)              #2*t
    twot_inv = crad_inv(twot_inv)       #1/(2*t)
    s = fp_add(tmp, tmp)                #2*A + 2*r
    s = fp_add(s, r)                    #2*A + 3*r
    s = fp_mul(r,s)                     #r*(2*A + 3*r)
    s = fp_add(1,s)                     #1 + r*(2*A + 3*r)
    s = fp_mul(s, twot_inv)             #(1 + r*(2*A + 3*r))/(2*t)
    
    u_inv = 1

    M = fp_mul(s, u_inv)                #s*u_inv
    M = fp_add(M,M)                     #2*s*u_inv
    N = fp_mul(u_inv, u_inv)
    N = fp_mul(u_inv, N)                #u_inv^3
    N = fp_mul(t, N)                    #t*u_inv^3
    N = fp_add(N,N)                     #2*t*u_inv^3
    
    return M, N

def three_iso(A, B):
    '''
    Adapted from the Radical isogenies code: 
    Applies a single radical isogeny of degree 3 by mapping A, B to newA, newB.
    A, B represent the curve y^2 + A*x*y + B*y = x^3,
    newA and newB represent the curve y^2 + newA*x*y + newB*y = x^3 
    (obtained from the 3-isogeny with  kernel <(0,0)>). 
    Due to the limitations of constant-time implementations, we obtain no speed-up
    with 9-isogenies, so we use 3-isogenies. As there is no inversion, we do not need
    a projective version.
    '''
    C = crad_tri(B)
    AC = fp_mul(A,C)
    
    tmp = fp_add(C,C)
    tmp = fp_add(C, tmp)            #3C
    newA = fp_add(tmp, tmp)         #6C
    newA = fp_sub(A, newA)          #A - 6C
    
    tmpb = fp_add(B,B)
    tmpb = fp_add(tmpb, tmpb)
    tmpb = fp_add(tmpb, tmpb)       #8B
    tmpb = fp_add(B, tmpb)          #9B
    
    newB = fp_sub(A, tmp)           #A - 3C
    newB = fp_mul(AC, newB)         #A*AC - 3*AC*C
    newB = fp_add(newB, tmpb)       #A*AC - 3*AC*C + 9B

    return newA, newB

def act_with_three_on_Montgomery(A, exp, Tp_proj = None):
    '''
    Adapted from the Radical isogenies code.
    Applies exp number of 3 isogenies to the Montgomery curve E_A on the floor.
    The sign of exp indicates the direction of the isogenies.
    Due to the limitations of constant-time implementations, we obtain no speed-up
    with 9-isogenies, so we use 3-isogenies. As there is no inversion, we do not need
    a projective version.
    '''
    
    A = (A * sign(exp)) % p         #
    if Tp_proj == None:
        A_proj = [fp_add(A, 2), 4]
        Tp_proj, _ = sampling_ell_order_point(A, 0)
        Tp_proj_times_3 = xMUL(Tp_proj, A_proj, 0)
        while isinfinity(Tp_proj) and isinfinity(Tp_proj_times_3):
            Tp_proj, _ = sampling_ell_order_point(A, 0)
            Tp_proj_times_3 = xMUL(Tp_proj, A_proj, 0)

        if not isinfinity(Tp_proj_times_3):
            Tp_proj = list(Tp_proj_times_3)

        sign_exp = 1
    else:
        sign_exp = sign(exp)

    # Affine x-coordinate Tp
    xTp = crad_inv(Tp_proj[1])
    xTp = fp_mul(xTp, Tp_proj[0])
    xTp = (xTp * sign_exp) % p         # <--- isomorphism from A to -A with zera-trace points on A
    
    M, N = Montgomery_to_Tate_three(A, xTp)         #<--- confusing line of code in magma
    M = (-M) % p
    
    for i in range(0, abs(exp)):
        M, N = three_iso(M,N)

    for i in range(abs(exp), e_3, 1):        # Exponent bound of 3 is assumed to be e_3, thus c_3 degree-3 radical isogenies
        B0, B1 = three_iso(M,N)
    
    M = (-M) % p
    
    A = Weier_to_Montgomery([M,0,N,0,0])
    A = (A * sign(exp)) % p         #
    
    return A

def nine_iso(A):
    '''
    Adapted from the Radical isogenies code: 
    Applies a single radical isogeny of degree 9 by mapping A to newA (output).
    Due to the limitations of constant-time implementations, we obtain no speed-up
    with 9-isogenies, so this (affine) radical isogeny is not used. 
    '''


    A2 = fp_sqr(A)              #A^2
    A3 = fp_mul(A2, A)          #A^3
    A4 = fp_sqr(A2)             #A^4
    A5 = fp_mul(A4, A)          #A^5
    A6 = fp_sqr(A3)             #A^6
    A7 = fp_mul(A6, A)          #A^7
    A8 = fp_sqr(A4)             #A^8
    A1_2 = fp_add(A,A)          #2*A
    A1_3 = fp_add(A1_2,A)       #3*A etc etc
    A1_4 = fp_add(A1_2,A1_2)
    A1_5 = fp_add(A1_2,A1_3)
    A1_6 = fp_add(A1_3,A1_3)
    A1_10 = fp_add(A1_6,A1_4)
    A1_23 = fp_mul(23,A)
    A2_2 = fp_add(A2,A2)
    A2_3 = fp_add(A2,A2_2)
    A2_4 = fp_add(A2_2,A2_2)
    A2_5 = fp_add(A2_3,A2_2)
    A2_8 = fp_add(A2_4,A2_4)
    A2_12 = fp_add(A2_8,A2_4)
    A2_16 = fp_add(A2_8,A2_8)
    A2_56 = fp_mul(7,A2_8)
    A3_2 = fp_add(A3,A3)
    A3_3 = fp_add(A3_2,A3)
    A3_6 = fp_add(A3_3,A3_3)
    A3_8 = fp_add(A3_6,A3_2)
    A3_10 = fp_add(A3_8,A3_2)
    A3_14 = fp_add(A3_8,A3_6)
    A3_16 = fp_add(A3_8,A3_8)
    A3_19 = fp_add(A3_16,A3_3)
    A3_89 = fp_mul(89, A3)
    A4_2 = fp_add(A4, A4)
    A4_3 = fp_add(A4_2,A4)
    A4_5 = fp_add(A4_2,A4_3)
    A4_6 = fp_add(A4_3,A4_3)
    A4_8 = fp_add(A4_2,A4_6)
    A4_12 = fp_add(A4_6,A4_6)
    A4_16 = fp_add(A4_8,A4_8)
    A4_20 = fp_add(A4_12,A4_8)
    A4_98 = fp_mul(98, A4)
    A5_2 = fp_add(A5, A5)
    A5_3 = fp_add(A5, A5_2)
    A5_6 = fp_add(A5_3, A5_3)
    A5_15 = fp_mul(15, A5)
    A5_77 = fp_mul(77, A5)
    A6_2 = fp_add(A6, A6)
    A6_7 = fp_mul(7, A6)
    A6_41 = fp_mul(41, A6)
    A7_2 = fp_add(A7, A7)
    A7_14 = fp_mul(14, A7)
    A8_2 = fp_add(A8, A8)

    C = fp_sub(A5, A4)                  #A5 - A4
    C_tmp = fp_sub(A2, A)               #A2 - A
    C_tmp = fp_add(C_tmp, 1)            #A2 - A + 1
    C = fp_mul(C, C_tmp)                #(A5 - A4)*(A2 - A + 1)
    C = fp_mul(C, C_tmp)                #(A5 - A4)*(A2 - A + 1)^2
    C = fp_mul(C, C_tmp)                #(A5 - A4)*(A2 - A + 1)^3
    C = crad_novem(C)                   #((A5 - A4)*(A2 - A + 1)^3)^novem_exp

    inv_tmp = fp_mul(C_tmp, A)          #(A3 - A2 + A)
    inv = fp_sqr(inv_tmp)               #(A3 - A2 + A)^2
    inv = fp_mul(inv, inv_tmp)          #(A3 - A2 + A)^3
    inv_tmp = fp_sub(inv_tmp, A2_5)     #(A3 - 6*A2 + A)
    inv_tmp = fp_add(inv_tmp, A1_2)     #(A3 - 6*A2 + 3*A)
    inv_tmp = fp_add(inv_tmp, 1)        #(A3 - 6*A2 + 3*A + 1)
    inv = fp_mul(inv, inv_tmp)          #(A3 - A2 + A)^3*(A3 - 6*A2 + 3*A + 1)
    inv = crad_inv(inv)                 #1/((A3 - A2 + A)^3*(A3 - 6*A2 + 3*A + 1))

    S1 = fp_sub(A2_2, A1_5)
    S1 = fp_add(S1, 2)                  #S1 = (2*A2 - 5*A + 2)

    S2 = fp_add(A3, 1)
    S2 = fp_mul(-1,S2)                  #S2 = (-A3 - 1)

    S3 = fp_sub(A1_3, A2_3)
    S3 = fp_sub(S3, 3)
    S3 = fp_mul(S3, A2)                 #S3 = (-3*A2 + 3*A - 3)*A2

    S4 = fp_sub(A5, A4_6)
    S4 = fp_add(S4, A3_8)
    S4 = fp_sub(S4, A2_8)
    S4 = fp_add(S4, A1_3)
    S4 = fp_sub(S4, 1)
    S4 = fp_mul(S4, A)                  #S4 = (A5 - 6*A4 + 8*A3 - 8*A2 + 3*A - 1)*A

    S5 = fp_sub(A5_2, A4_8)
    S5 = fp_add(S5, A3_14)
    S5 = fp_sub(S5, A2_16)
    S5 = fp_add(S5, A1_10)
    S5 = fp_sub(S5, 4)
    S5 = fp_mul(S5, A2)                 #S5 = (2*A5 - 8*A4 + 14*A3 - 16*A2 + 10*A - 4)*A2

    S6 = fp_sub(A6, A5_6)
    S6 = fp_add(S6, A4_12)
    S6 = fp_sub(S6, A3_16)
    S6 = fp_add(S6, A2_12)
    S6 = fp_sub(S6, A1_6)
    S6 = fp_add(S6, 1)
    S6 = fp_mul(S6, A2)                 #S6 = (A6 - 6*A5 + 12*A4 - 16*A3 + 12*A2 - 6*A + 1)*A2

    S7 = fp_sub(A5_3, A6_2)
    S7 = fp_sub(S7, A4_3)
    S7 = fp_sub(S7, A3)
    S7 = fp_add(S7, A2_3)
    S7 = fp_sub(S7, A1_3)
    S7 = fp_add(S7, 1)
    S7 = fp_mul(S7, A3)                 #S7 = (-2*A6 + 3*A5 - 3*A4 - A3 + 3*A2 - 3*A + 1)*A3

    S8 = fp_sub(A6_7, A7_2)
    S8 = fp_sub(S8, A5_15)
    S8 = fp_add(S8, A4_20)
    S8 = fp_sub(S8, A3_19)
    S8 = fp_add(S8, A2_12)
    S8 = fp_sub(S8, A1_5)
    S8 = fp_add(S8, 1)
    S8 = fp_mul(S8, A3)                 #S8 = (-2*A7 + 7*A6 - 15*A5 + 20*A4 - 19*A3 + 12*A2 - 5*A + 1)*A3

    S9 = fp_sub(A8_2, A7_14)
    S9 = fp_add(S9, A6_41)
    S9 = fp_sub(S9, A5_77)
    S9 = fp_add(S9, A4_98)
    S9 = fp_sub(S9, A3_89)
    S9 = fp_add(S9, A2_56)
    S9 = fp_sub(S9, A1_23)
    S9 = fp_add(S9, 5)
    S9 = fp_mul(S9, A4)                 #S9 = (2*A8 - 14*A7 + 41*A6 - 77*A5 + 98*A4 - 89*A3 + 56*A2 - 23*A + 5)*A4

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
    output = fp_mul(output, C)
    output = fp_add(output, S7)
    output = fp_mul(output, C)
    output = fp_add(output, S8)
    output = fp_mul(output, C)
    output = fp_add(output, S9)

    output = fp_mul(inv, output)        #see crad magma for explanation on these lines

    return output

def nine_iso_projective(X, Z):
    '''
    Adapted from the Radical isogenies code: 
    Applies a single radical isogeny of degree 9 by mapping A to newA (output).
    Due to the limitations of constant-time implementations, we obtain no speed-up
    with 9-isogenies, so this (affine) radical isogeny is not used. 
    '''


    X2 = fp_sqr(X)              #X^2
    X3 = fp_mul(X2, X)          #X^3
    X4 = fp_sqr(X2)             #X^4
    X5 = fp_mul(X4, X)          #X^5
    X6 = fp_sqr(X3)             #X^6
    X7 = fp_mul(X6, X)          #X^7
    X8 = fp_sqr(X4)             #X^8
    
    Z2 = fp_sqr(Z)              #Z^2
    Z3 = fp_mul(Z2, Z)          #Z^3
    Z4 = fp_sqr(Z2)             #Z^4
    Z5 = fp_mul(Z4, Z)          #Z^5
    Z6 = fp_sqr(Z3)             #Z^6
    Z7 = fp_mul(Z6, Z)          #Z^7
    Z8 = fp_sqr(Z4)             #Z^8
    
    XZ = fp_mul(X, Z)
    
    X2Z = fp_mul(X2, Z)
    XZ2 = fp_mul(X, Z2)
    
    X4Z = fp_mul(X4, Z)
    X3Z2 = fp_mul(X3, Z2)
    X2Z3 = fp_mul(X2, Z3)
    XZ4 = fp_mul(X, Z4)
    
    X5Z = fp_mul(X5, Z)
    X4Z2 = fp_mul(X4, Z2)
    X3Z3 = fp_mul(X3, Z3)
    X2Z4 = fp_mul(X2, Z4)
    XZ5 = fp_mul(X, Z5)
    
    X6Z = fp_mul(X6, Z)
    X5Z2 = fp_mul(X5, Z2)
    X4Z3 = fp_mul(X4, Z3)
    X3Z4 = fp_mul(X3, Z4)
    X2Z5 = fp_mul(X2, Z5)
    XZ6 = fp_mul(X, Z6)
    
    X7Z = fp_mul(X7, Z)
    X6Z2 = fp_mul(X6, Z2)
    X5Z3 = fp_mul(X5, Z3)
    X4Z4 = fp_mul(X4, Z4)
    X3Z5 = fp_mul(X3, Z5)
    X2Z6 = fp_mul(X2, Z6)
    XZ7 = fp_mul(X, Z7)
     
    T1 = fp_mul(X4, Z)
    T1 = fp_sub(X5, T1)                  #X5 - X4*Z
    T2 = fp_sub(X2, XZ)
    T2 = fp_add(T2, Z2)                 #X2 - XZ + Z2
    T23 = fp_sqr(T2)
    T23 = fp_mul(T23, T2)
    
    C = fp_mul(T1, T23)
    C = fp_mul(Z7, C)                   #Z^7 * (X^5 - X^4*Z) (X^2 - X*Z + Z^2)^3
    C = fp_exp(C, novem_exp)
    
    inv_T1 = fp_sub(X3, X2Z)
    inv_T1 = fp_add(inv_T1, XZ2)
    inv_T13 = fp_sqr(inv_T1)
    inv_T13 = fp_mul(inv_T13, inv_T1)       #(X^3 - X^2*Z + X*Z^2)^3
    
    inv_T2 = fp_sub(X3, fp_mul(6, X2Z))
    inv_T2 = fp_add(inv_T2, fp_mul(3, XZ2))
    inv_T2 = fp_add(inv_T2, Z3)             #(X^3 - 6x^2*z + 3XZ^2 + Z^3) 
    
    inv = fp_mul(inv_T13, inv_T2)
    inv = fp_mul(Z6, inv)

    #S1 = S1p/Z^2
    S1p = fp_add(X2, X2)
    S1p = fp_sub(S1p, fp_mul(5, XZ))
    S1p = fp_add(S1p, fp_add(Z2, Z2))       #(2X^2 - 5XZ + 2Z)
    
    #S2 = S2p/Z^3
    S2p = fp_add(X3, Z3)
    S2p = fp_mul(-1, S2p)                   #(-X^3 - Z^3)
    
    #S3 = S3p/Z^4
    S3p = fp_sub(XZ, X2)
    S3p = fp_sub(S3p, Z2)
    S3p = fp_mul(S3p, X2)
    S3p = fp_add(S3p, fp_add(S3p, S3p))     #(-X^2 + XZ - Z^2)*3X^2
    
    #S4 = S4p/Z^6    
    S4p = fp_sub(X5, fp_mul(6, X4Z))
    S4p = fp_add(S4p, fp_mul(8, X3Z2))
    S4p = fp_add(S4p, fp_mul(-8, X2Z3))
    S4p = fp_add(S4p, fp_mul(3, XZ4))
    S4p = fp_sub(S4p, Z5)
    S4p = fp_mul(S4p, X)                    #(X5 - 6*X^4*Z + 8*X^3*Z^2 - 8*X^2*Z^3 + 3*X*Z^4 - Z^5)*X
    
    #S5 = S5p/Z^7
    S5p = fp_sub(fp_add(X5, X5), fp_mul(8, X4Z))
    S5p = fp_add(S5p, fp_mul(14, X3Z2))
    S5p = fp_add(S5p, fp_mul(-16, X2Z3))
    S5p = fp_add(S5p, fp_mul(10, XZ4))
    S5p = fp_sub(S5p, fp_add(fp_add(Z5, Z5), fp_add(Z5, Z5)))
    S5p = fp_mul(S5p, X2)
    
    #S6 = S6p/Z^8
    S6p = fp_sub(X6, fp_mul(6, X5Z))
    S6p = fp_add(S6p, fp_mul(12, X4Z2))
    S6p = fp_add(S6p, fp_mul(-16, X3Z3))
    S6p = fp_add(S6p, fp_mul(12, X2Z4))
    S6p = fp_add(S6p, fp_mul(-6, XZ5))
    S6p = fp_add(S6p, Z6)
    S6p = fp_mul(S6p, X2)                #(X6 - 6*X5Z + 12*X4Z2 - 16*X3Z3 + 12*X2Z4 - 6*XZ5 + Z6)*X2
    
    #S7 = S7p/Z^9
    S7p = fp_add(X6, X6)
    S7p = fp_sub(fp_mul(3, X5Z), S7p)
    S7p = fp_add(S7p, fp_mul(-3, X4Z2))
    S7p = fp_sub(S7p, X3Z3)
    S7p = fp_add(S7p, fp_mul(3, X2Z4))
    S7p = fp_add(S7p, fp_mul(-3, XZ5))
    S7p = fp_add(S7p, Z6)
    S7p = fp_mul(S7p, X3)               #(-2*X6 + 3*X5Z - 3*X4Z2 - X3Z3 + 3*X2Z4 - 3*XZ5 + Z6)*X3

    #S8 = S8p/Z^10
    S8p = fp_add(X7, X7)
    S8p = fp_sub(fp_mul(7, X6Z), S8p)
    S8p = fp_add(S8p, fp_mul(-15, X5Z2))
    S8p = fp_add(S8p, fp_mul( 20, X4Z3))
    S8p = fp_add(S8p, fp_mul(-19, X3Z4))
    S8p = fp_add(S8p, fp_mul( 12, X2Z5))
    S8p = fp_add(S8p, fp_mul(-5,  XZ6))
    S8p = fp_add(S8p, Z7)
    S8p = fp_mul(S8p, X3)               #(-2*X7 + 7*X6Z - 15*X5Z2 + 20*X4Z3 - 19*X3Z4 + 12*X2Z5 - 5*XZ6 + Z7)*X3
    
    #S9 = S9p/Z^12
    S9p = fp_add(X8, X8)
    S9p = fp_add(S9p, fp_mul(-14, X7Z))
    S9p = fp_add(S9p, fp_mul( 41, X6Z2))
    S9p = fp_add(S9p, fp_mul(-77, X5Z3))
    S9p = fp_add(S9p, fp_mul( 98, X4Z4))
    S9p = fp_add(S9p, fp_mul(-89, X3Z5))
    S9p = fp_add(S9p, fp_mul( 56,  X2Z6))
    S9p = fp_add(S9p, fp_mul(-23, XZ7))
    S9p = fp_add(S9p, fp_mul(  5, Z8))
    S9p = fp_mul(S9p, X4)               #(2*X8 - 14*X7Z + 41*X6Z2 - 77*X5Z3 + 98*X4Z4 - 89*X3Z5 + 56*X26 - 23*XZ7 + 5Z8)*X4


    output = fp_mul(S1p, C)
    output = fp_add(output, fp_mul(Z, S2p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z2, S3p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z2, S4p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z3, S5p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z4, S6p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z5, S7p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z6, S8p))
    output = fp_mul(output, C)
    output = fp_add(output, fp_mul(Z6, S9p))
    
    return output, inv

def pro_to_aff_nine(X,Z):
    A = fp_inv(Z)
    A = fp_mul(X, A)
    
    return A

def act_with_nine_on_Montgomery(A, exp, Tp_proj = None):
    '''
    Adapted from the Radical isogenies code.
    Applies exp number of 9 isogenies to the Montgomery curve E_A on the floor.
    The sign of exp indicates the direction of the isogenies.
    Due to the limitations of constant-time implementations, we obtain no speed-up
    with 9-isogenies, so this function is not used. 
    '''
    
    A = (A * sign(exp)) % p         #
    if Tp_proj == None:
        Tp_proj, _ = sampling_ell_order_point(A, 0)
        while isinfinity(xMUL(Tp_proj, [fp_add(A, 2), 4], 0)):
            Tp_proj, _ = sampling_ell_order_point(A, 0)

        sign_exp = 1
    else:
        sign_exp = sign(exp)

    # Affine x-coordinate Tp
    xTp = crad_inv(Tp_proj[1])
    xTp = fp_mul(xTp, Tp_proj[0])
    xTp = (xTp * sign_exp) % p         # <--- isomorphism from A to -A with zera-trace points on A
    
    M, N = Montgomery_to_Tate(A, xTp)

    A_inv = fp_sub(M, N)
    A_inv = fp_sub(A_inv, 1)        
    A_inv = crad_inv(A_inv)           #1/(M-N-1)

    A = fp_sub(M, 1)
    A = fp_sqr(A)                   #(M - 1)^2
    A = fp_mul(A, A_inv)            #(M - 1)^2/(M-N-1)
    
    for i in range(0, abs(exp), 1):         #no -1 because no last velu
        A = nine_iso(A)

    for i in range(abs(exp), e_3 // 2, 1):        # Exponent bound of 3 is assumed to be e_3, thus c_3 // 2 degree-9 radical isogenies
        dummy = nine_iso(A)

    A2 = fp_sqr(A)
    A3 = fp_mul(A2,A)
    A4 = fp_sqr(A2)
    A5 = fp_mul(A4,A)
    A6 = fp_sqr(A3)
    A7 = fp_mul(A6,A)
    A8 = fp_sqr(A4)
    A9 = fp_mul(A8,A)
    A10 = fp_sqr(A5)
    A11 = fp_mul(A10,A)
    A12 = fp_sqr(A6)
    A13 = fp_mul(A12,A)
    A14 = fp_sqr(A7)
    A15 = fp_mul(A14,A)
    A16 = fp_sqr(A8)
    A17 = fp_mul(A16,A)

    a1 = fp_sub(A2, A3)
    a1 = fp_add(a1, 1)                     #-A3 + A2 + 1

    a2 = fp_sub(A4, A3)
    a2 = fp_add(a2, a2)
    a2 = fp_sub(a2, A5)
    a2 = fp_add(a2, A2)                    #-A5 + 2*A4 - 2*A3 + A2
    
    a3 = a2

    a4 = 0                                  #we apply no last Velu
    
    a6 = 0                                  #we apply no last Velu
    
    # a4 = fp_sub(A10, A11)
    # a4 = fp_add(a4, fp_mul(   8, A9))
    # a4 = fp_add(a4, fp_mul( -33, A8))
    # a4 = fp_add(a4, fp_mul(  72, A7))
    # a4 = fp_add(a4, fp_mul(-108, A6))
    # a4 = fp_add(a4, fp_mul( 114, A5))
    # a4 = fp_add(a4, fp_mul( -81, A4))
    # a4 = fp_add(a4, fp_mul(  37, A3))
    # a4 = fp_add(a4, fp_mul( -10, A2))
    # a4 = fp_add(a4, A)
    # a4 = fp_mul(5, a4)                      #(-5*A11 + 5*A10 + 40*A9 - 165*A8 + 360*A7 - 540*A6 + 570*A5 - 405*A4 + 185*A3 - 50*A2 + 5*A)

    # a6 = fp_sub(A, A17)
    # a6 = fp_add(a6, fp_mul(    -7, A16))
    # a6 = fp_add(a6, fp_mul(    63, A15))
    # a6 = fp_add(a6, fp_mul(  -230, A14))
    # a6 = fp_add(a6, fp_mul(   641, A13))
    # a6 = fp_add(a6, fp_mul( -1639, A12))
    # a6 = fp_add(a6, fp_mul(  3691, A11))
    # a6 = fp_add(a6, fp_mul( -6707, A10))        
    # a6 = fp_add(a6, fp_mul(  9425,  A9))
    # a6 = fp_add(a6, fp_mul(-10174,  A8))
    # a6 = fp_add(a6, fp_mul(  8456,  A7))
    # a6 = fp_add(a6, fp_mul( -5379,  A6))
    # a6 = fp_add(a6, fp_mul(  2559,  A5))
    # a6 = fp_add(a6, fp_mul(  -865,  A4))
    # a6 = fp_add(a6, fp_mul(   190,  A3))
    # a6 = fp_add(a6, fp_mul(   -24,  A2))    # (-A17 - 7*A16 + 63*A15 - 230*A14 + 641*A13 - 1639*A12 + 3691*A11 - 6707*A10 + 9425*A9 - 10174*A8 + 8456*A7 - 5379*A6 + 2559*A5 - 865*A4 + 190*A3 - 24*A2 + A);

    output = Weier_to_Montgomery([a1,a2,a3,a4,a6])
    output = (output * sign(exp)) % p

    return output