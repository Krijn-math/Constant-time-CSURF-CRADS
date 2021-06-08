#from src.radical.crad_general import *
#from src.montgomery import isinfinity

'''material needed for radical 3 and 9-isogenies'''
'''computations verified against Magma implementation CRAD'''

#e_3 = 16

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
    t0 = crad_sq(t0)
   
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



def Montgomery_to_Tate_nine_pro_with_affine_point(A, xP):
    '''
    This is a projective version of the original Montgomery to Tate function
    specifically adapted to provide the right projective coordinates for the
    degree-nine projective radical isogeny. This assumes the point is given
    in affine coordinates.
'''

    r = xP
    tmp = fp_add(A,r)                   #A+r
    t = fp_mul(r,tmp)                   #r*(A+r)
    t = fp_add(1, t)                    #(1 + r*(A+r))
    t = fp_mul(r,t)                     #r * (1 + r*(A+r))
    t = crad_sq(t)                      #(r*(1+r*(A+r)))^sq_exp
    
    twot = fp_add(t,t)                  #2*t
    twot2 = fp_sqr(twot)                #(2*t)^2
    twot4 = fp_sqr(twot2)                #(2*t)^4
    twot8 = fp_sqr(twot4)                #(2*t)^8

    s = fp_add(tmp, tmp)                #2*A + 2*r
    s = fp_add(s, r)                    #2*A + 3*r
    s = fp_mul(r,s)                     #r*(2*A + 3*r)
    s = fp_add(1,s)                     #1 + r*(2*A + 3*r)
    
    u = fp_add(r,r)                     #2*r
    u = fp_add(tmp, u)                  #A + 3*r
    u = fp_mul(u, twot2)                #(2t)^2(A + 3*r)
    s_sq = fp_sqr(s)                    #s^2
    u = fp_sub(u, s_sq)                 #(2t)^2*A + (2t)^2*3*r - s^2
    
    u3 = fp_sqr(u)
    u3 = fp_mul(u3, u)
    
    Mt = fp_mul(s, u)
    Mt = fp_add(Mt, Mt)                 #2*s*u
    
    M = fp_sub(Mt, twot4)
    M = fp_sqr(M)                       #(Mt - (2t)^4)^2
    
    N = fp_mul(twot4, Mt)
    N = fp_sub(N, u3)
    N = fp_sub(N, twot8)                #(2t)^4 - u^3 - (2t)^8
    
    return M, N

def Weier_to_Montgomery_pro(coeffs):
    '''
   This is a projective version of the Weier to Montgomery function, but
   assumes for simplicity that coeffs 3 and 4 are 0.
   '''
    a1p = coeffs[0]
    a2p = coeffs[1]
    a3p = coeffs[2]
    a4p = coeffs[3]
    a6p = coeffs[4]
    X = coeffs[5]
    Z = coeffs[6]
    Z2 = fp_sqr(Z)
    Z4 = fp_sqr(Z2)
    Z6 = fp_mul(Z4, Z2)
    Z8 = fp_sqr(Z4)

    sp = fp_sub(0,a1p)               #-a1p

    Bp = fp_sqr(a1p)
    Bp = fp_add(Bp, fp_mul(a2p, Z))
    Bp = fp_add(Bp, fp_mul(a2p, Z))
    Bp = fp_add(Bp, fp_mul(a2p, Z))
    Bp = fp_add(Bp, fp_mul(a2p, Z))     #a1p^2 + 4*a2p*Z

    Cp = fp_mul(a1p, a3p)               #a1p*a3p

    Dp = fp_sqr(a3p)                      #a3p^2
    
    d_tmp = fp_mul(Cp, Z4)
    d_tmp = fp_add(d_tmp, d_tmp)                    #2*Cp*Z^4             
    d_tmp = fp_add(d_tmp, d_tmp)                    #4*Cp*Z^4  
    d_tmp = fp_add(d_tmp, d_tmp)                    #8*Cp*Z^4
    d_tmp = fp_add(d_tmp, fp_add(d_tmp, d_tmp))     #24*Cp*Z^4

    Bp2 = fp_sqr(Bp)
    d0p = fp_sub(Bp2, d_tmp)                        #d0p = Bp^2 - 24*Cp*Z^4 

    d_tmp = fp_add(d_tmp, fp_add(d_tmp, d_tmp))     #72*Cp*Z^4  
    
    d2_tmp = fp_mul(Z8, Dp)
    d2_tmp = fp_add(d2_tmp, d2_tmp)                    #2*Dp*Z^8             
    d2_tmp = fp_add(d2_tmp, d2_tmp)                    #4*Dp*Z^8  
    d2_tmp = fp_add(d2_tmp, d2_tmp)                    #8*Dp*Z^8  
    d2_tmp = fp_add(d2_tmp, fp_add(d2_tmp, d2_tmp))     #24*Dp*Z^8   
    d2_tmp = fp_add(d2_tmp, fp_add(d2_tmp, d2_tmp))     #72*Dp*Z^8
    d2_tmp = fp_add(d2_tmp, fp_add(d2_tmp, d2_tmp))     #72*Dp*Z^8
    d2_tmp = fp_add(d2_tmp, d2_tmp)                     #4^2*3^3*Z^8*Dp

    d1p = fp_add(Bp2, Bp2)
    d1p = fp_sub(d1p, d_tmp)
    d1p = fp_mul(Bp, d1p)
    d1p = fp_add(d1p, d2_tmp)                           #d1p = Bp(2*Bp^2 - 72*Z^4*Cp) + 4^2*3^3*Z^8*Dp

    d0p3 = fp_sqr(d0p)
    d0p3 = fp_mul(d0p, d0p3)
    rootp = fp_sqr(d1p)
    rootp = fp_sub(rootp, d0p3)
    rootp = fp_sub(rootp, d0p3)
    rootp = fp_sub(rootp, d0p3)
    rootp = fp_sub(rootp, d0p3)
    rootp = crad_sq(rootp)                   #(d1p^2-4*d0p^3)^sq_exp

    cbp = fp_add(d1p, rootp)
    cbp = fp_mul(cbp, inv_2)
    cbp = crad_tri(cbp)                   #((d1p+rootp)/2)^tri_exp;

    rp = fp_add(Bp, cbp)
    rp = fp_mul(-cbp, rp)
    rp = fp_sub(rp, d0p)                        #rp = -cbp(Bp + cbp) - d0p
    
    rb = fp_mul(Z6, cbp)
    rb = fp_add(rb, rb)
    rb = fp_add(rb, rb)
    rb = fp_add(rb, fp_add(rb,rb))          #3*4*Z^6*cbp
    
    t_tmp = fp_mul(cbp, a3p)             
    t_tmp = fp_mul(t_tmp, Z4)
    t_tmp = fp_add(t_tmp, t_tmp)
    t_tmp = fp_add(t_tmp, t_tmp)
    t_tmp = fp_add(t_tmp, fp_add(t_tmp, t_tmp))     #3*4*Z^4*cbp*a3p
    
    tp = fp_mul(rp, a1p)
    tp = fp_add(t_tmp, tp)
    tp = fp_sub(0, tp)                              # -(3*4*cbp*Z^4*a3p + rp*a1p)
     
    _2cbp = fp_add(cbp, cbp)
    _3cbp = fp_add(cbp, _2cbp)
    _6cbp = fp_add(_3cbp, _3cbp)
    _12cbp = fp_add(_6cbp, _6cbp)
    _24cbp = fp_add(_12cbp, _12cbp)
    _72cbp = fp_add(_24cbp, fp_add(_24cbp, _24cbp))
       
    AZ1 = fp_mul(_72cbp, cbp)
    AZ1 = fp_mul(AZ1, Z4)
    AZ1 = fp_mul(AZ1, sp)
    AZ1 = fp_mul(AZ1, a3p)                  #72*Z^4*cbp^2*sp*a3p
    
    AZ2 = fp_mul(Z, _24cbp)
    AZ2 = fp_mul(AZ2, rp)
    AZ2 = fp_mul(AZ2, a2p)                  #3*8*Z*cbp*(rp*a2p)

    AZ3 = fp_mul(rp, sp)
    AZ3 = fp_add(tp, AZ3)
    AZ3 = fp_mul(AZ3, a1p)
    AZ3 = fp_mul(AZ3, _6cbp)                #(tp + rp*sp)a1p*(3*2*cbp)

    
    AZ4 = fp_sqr(rp)
    AZ4 = fp_add(AZ4, fp_add(AZ4, AZ4))     #3*rp2
    
    AZ5 = fp_mul(_6cbp, fp_mul(sp, tp))     #6*cbp*sp*tp

    AZ = fp_sub(AZ2, AZ1)
    AZ = fp_sub(AZ, AZ3)
    AZ = fp_add(AZ, AZ4)
    AZ = fp_sub(AZ, AZ5)
    AZ = fp_sqr(crad_quart(fp_mul(AZ, fp_sqr(rb))))      #in order to get correct root
    
    AX_tmp = fp_add(a1p, a1p)
    AX_tmp = fp_add(AX_tmp, sp)
    AX_tmp = fp_mul(sp, AX_tmp)

    AX = fp_mul(Z, a2p)
    AX = fp_add(AX, AX)
    AX = fp_add(AX, AX)
    AX = fp_sub(AX, AX_tmp)
    AX = fp_mul(cbp, AX)
    AX = fp_add(rp, AX)
    AX2 = fp_add(AX, AX)
    AX = fp_add(AX2, AX)
    AX = fp_mul(AX, rb)                                 #in order to get correct root
    
    return AX, AZ

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
    C = crad_novem(C)
    
    inv_T1 = fp_sub(X3, X2Z)
    inv_T1 = fp_add(inv_T1, XZ2)
    inv_T13 = fp_sqr(inv_T1)
    inv_T13 = fp_mul(inv_T13, inv_T1)       #(X^3 - X^2*Z + X*Z^2)^3
       
    _2X2Z = fp_add(X2Z, X2Z)
    _6X2Z = fp_add(_2X2Z, fp_add(_2X2Z, _2X2Z))
    _3XZ2 = fp_add(XZ2, fp_add(XZ2, XZ2))
    
    inv_T2 = fp_sub(X3, _6X2Z)
    inv_T2 = fp_add(inv_T2, _3XZ2)
    inv_T2 = fp_add(inv_T2, Z3)             #(X^3 - 6x^2*z + 3XZ^2 + Z^3) 
    
    inv = fp_mul(inv_T13, inv_T2)
    inv = fp_mul(Z6, inv)
    
    _2XZ = fp_add(XZ, XZ)
    _4XZ = fp_add(_2XZ, _2XZ)
    _5XZ = fp_add(_4XZ, XZ)
    
    #S1 = S1p/Z^2
    S1p = fp_add(X2, X2)
    S1p = fp_sub(S1p, _5XZ)
    S1p = fp_add(S1p, fp_add(Z2, Z2))       #(2X^2 - 5XZ + 2Z)
    
    #S2 = S2p/Z^3
    S2p = fp_add(X3, Z3)
    S2p = fp_sub(0, S2p)                   #(-X^3 - Z^3)
    
    #S3 = S3p/Z^4
    S3p = fp_sub(XZ, X2)
    S3p = fp_sub(S3p, Z2)
    S3p = fp_mul(S3p, X2)
    S3p = fp_add(S3p, fp_add(S3p, S3p))     #(-X^2 + XZ - Z^2)*3X^2
    
    #S4 = S4p/Z^6    
    _2X4Z = fp_add(X4Z, X4Z)
    _6X4Z = fp_add(_2X4Z, fp_add(_2X4Z, _2X4Z))
    _2X3Z2 = fp_add(X3Z2, X3Z2)
    _8X3Z2 = fp_add(fp_add(_2X3Z2, _2X3Z2), fp_add(_2X3Z2, _2X3Z2))
    _2X2Z3 = fp_add(X2Z3, X2Z3)
    _8X2Z3 = fp_add(fp_add(_2X2Z3, _2X2Z3), fp_add(_2X2Z3, _2X2Z3))
    _3XZ4 = fp_add(XZ4, fp_add(XZ4, XZ4))
    
    S4p = fp_sub(X5, _6X4Z)
    S4p = fp_add(S4p, _8X3Z2)
    S4p = fp_sub(S4p, _8X2Z3)
    S4p = fp_add(S4p, _3XZ4)
    S4p = fp_sub(S4p, Z5)
    S4p = fp_mul(S4p, X)                    #(X5 - 6*X^4*Z + 8*X^3*Z^2 - 8*X^2*Z^3 + 3*X*Z^4 - Z^5)*X
    
    #S5 = S5p/Z^7
    _8X4Z = fp_add(_2X4Z, _6X4Z)
    _14X3Z2 = fp_add(_8X3Z2, fp_sub(_8X3Z2, _2X3Z2))
    _16X2Z3 = fp_add(_8X2Z3, _8X2Z3)
    _4XZ4 = fp_add(_3XZ4, XZ4)
    _10XZ4 = fp_add(_4XZ4, fp_add(_3XZ4, _3XZ4))
    
    S5p = fp_sub(fp_add(X5, X5), _8X4Z)
    S5p = fp_add(S5p, _14X3Z2)
    S5p = fp_sub(S5p, _16X2Z3)
    S5p = fp_add(S5p, _10XZ4)
    S5p = fp_sub(S5p, fp_add(fp_add(Z5, Z5), fp_add(Z5, Z5)))
    S5p = fp_mul(S5p, X2)
    
    #S6 = S6p/Z^8
    _2X5Z = fp_add(X5Z, X5Z)
    _6X5Z = fp_add(_2X5Z, fp_add(_2X5Z, _2X5Z))
    _2X4Z2 = fp_add(X4Z2, X4Z2)
    _4X4Z2 = fp_add(_2X4Z2, _2X4Z2)
    _12X4Z2 = fp_add(_4X4Z2, fp_add(_4X4Z2, _4X4Z2))
    _2X3Z3 = fp_add(X3Z3, X3Z3)
    _8X3Z3 = fp_add(fp_add(_2X3Z3, _2X3Z3), fp_add(_2X3Z3, _2X3Z3))
    _16X3Z3 = fp_add(_8X3Z3, _8X3Z3)
    _2X2Z4 = fp_add(X2Z4, X2Z4)
    _4X2Z4 = fp_add(_2X2Z4, _2X2Z4)
    _12X2Z4 = fp_add(_4X2Z4, fp_add(_4X2Z4, _4X2Z4))
    _2XZ5 = fp_add(XZ5, XZ5)
    _6XZ5 = fp_add(_2XZ5, fp_add(_2XZ5, _2XZ5))
    
    S6p = fp_sub(X6, _6X5Z)
    S6p = fp_add(S6p, _12X4Z2)
    S6p = fp_sub(S6p, _16X3Z3)
    S6p = fp_add(S6p, _12X2Z4)
    S6p = fp_sub(S6p, _6XZ5)
    S6p = fp_add(S6p, Z6)
    S6p = fp_mul(S6p, X2)                #(X6 - 6*X5Z + 12*X4Z2 - 16*X3Z3 + 12*X2Z4 - 6*XZ5 + Z6)*X2
    
    #S7 = S7p/Z^9
    S7p = fp_add(X6, X6)
    S7p = fp_sub(fp_add(_2X5Z, X5Z), S7p)
    S7p = fp_sub(S7p, fp_add(_2X4Z2, X4Z2))
    S7p = fp_sub(S7p, X3Z3)
    S7p = fp_add(S7p, fp_add(_2X2Z4, X2Z4))
    S7p = fp_sub(S7p, fp_add(_2XZ5, XZ5))
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
    A = crad_inv(Z)
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
#    A_inv = crad_inv(A_inv)           #this step can be removed because of projective version

    A = fp_sub(M, 1)
    A = fp_sqr(A)                   #(M - 1)^2
#    A = fp_mul(A, A_inv)            #(M - 1)^2/(M-N-1)
    
    X, Z = A, A_inv
    
    for i in range(0, abs(exp), 1):         #no -1 because no last velu
        X, Z = nine_iso_projective(X,Z)

    for i in range(abs(exp), e_3 // 2, 1):        # Exponent bound of 3 is assumed to be e_3, thus c_3 // 2 degree-9 radical isogenies
        dummyX, dummyZ = nine_iso_projective(X,Z)

    A = pro_to_aff_nine(X, Z)

    A2 = fp_sqr(A)
    A3 = fp_mul(A2,A)
    A4 = fp_sqr(A2)
    A5 = fp_mul(A4,A)

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

def act_with_nine_on_Montgomery_pro_old(A, exp, Tp_proj = None):
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
    
    Tp_proj[0] = (Tp_proj[0] * sign_exp) % p
    
    #here, we use a specifically adapted version of Montgomery to Tate that
    #provides the right projective coordinates for the nine isogeny, which
    #saves the inversions used.
    X, Z = Montgomery_to_Tate_nine_pro(A, Tp_proj[0], Tp_proj[1]) 
    
    for i in range(0, abs(exp), 1):         #no -1 because no last velu
        X, Z = nine_iso_projective(X,Z)

    for i in range(abs(exp), e_3 // 2, 1):        # Exponent bound of 3 is assumed to be e_3, thus c_3 // 2 degree-9 radical isogenies
        dummyX, dummyZ = nine_iso_projective(X,Z)

    A = pro_to_aff_nine(X, Z)                   #TODO: we do not need to turn affine here, if we adjust Weier to Montg.

    A2 = fp_sqr(A)
    A3 = fp_mul(A2,A)
    A4 = fp_sqr(A2)
    A5 = fp_mul(A4,A)

    a1 = fp_sub(A2, A3)
    a1 = fp_add(a1, 1)                     #-A3 + A2 + 1
    
    a2 = fp_sub(A4, A3)
    a2 = fp_add(a2, a2)
    a2 = fp_sub(a2, A5)
    a2 = fp_add(a2, A2)                    #-A5 + 2*A4 - 2*A3 + A2
       
    a3 = a2

    a4 = 0                                  #we apply no last Velu
    
    a6 = 0                                  #we apply no last Velu
    

    output = Weier_to_Montgomery([a1,a2,a3,a4,a6])
    output = (output * sign(exp)) % p

    return output

def act_with_nine_on_Montgomery_pro(A, exp, Tp_proj = None):
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
    
    Tp_proj[0] = (Tp_proj[0] * sign_exp) % p
    
    #here, we use a specifically adapted version of Montgomery to Tate that
    #provides the right projective coordinates for the nine isogeny, which
    #saves the inversions used.
    X, Z = Montgomery_to_Tate_nine_pro(A, Tp_proj[0], Tp_proj[1]) 
    
    for i in range(0, abs(exp), 1):         #no -1 because no last velu
        X, Z = nine_iso_projective(X,Z)

    for i in range(abs(exp), e_3 // 2, 1):        # Exponent bound of 3 is assumed to be e_3, thus c_3 // 2 degree-9 radical isogenies
        dummyX, dummyZ = nine_iso_projective(X,Z)

    AX = X
    AZ = Z
    AX2 = fp_sqr(AX)
    AX3 = fp_mul(AX, AX2)
    AX4 = fp_sqr(AX2)
    AX5 = fp_mul(AX, AX4)
    
    AZ2 = fp_sqr(AZ)
    AZ3 = fp_mul(AZ, AZ2)
    AZ4 = fp_sqr(AZ2)
    
    a1p = fp_add(AZ2, AX2)
    a1p = fp_mul(Z, a1p)
    a1p = fp_sub(a1p, AX3)                          #-X^3 + X^2Z + Z^3
    

    
    AX4Z = fp_mul(AX4, AZ)
    AX3Z2 = fp_mul(AX3, AZ2)
    AX2Z3 = fp_mul(AX2, AZ3)
    
    a2p = fp_sub(AX2Z3, AX5)
    a2p = fp_add(a2p, fp_sub(AX4Z, AX3Z2))
    a2p = fp_add(a2p, fp_sub(AX4Z, AX3Z2))          #-X^5 + 2*X^4Z - 2*X^3Z^2 + X^2Z^3
  
    a3p = a2p

    # A = pro_to_aff_nine(X, Z)                   #TODO: we do not need to turn affine here, if we adjust Weier to Montg.

    # A2 = fp_sqr(A)
    # A3 = fp_mul(A2,A)
    # A4 = fp_sqr(A2)
    # A5 = fp_mul(A4,A)

    # a1 = fp_sub(A2, A3)
    # a1 = fp_add(a1, 1)                     #-A3 + A2 + 1

    # a2 = fp_sub(A4, A3)
    # a2 = fp_add(a2, a2)
    # a2 = fp_sub(a2, A5)
    # a2 = fp_add(a2, A2)                    #-A5 + 2*A4 - 2*A3 + A2
    
    # a3 = a2

    # a4 = 0                                  #we apply no last Velu
    
    # a6 = 0                                  #we apply no last Velu
    

    X, Z = Weier_to_Montgomery_pro([a1p,a2p,a3p,0,0, AX, AZ])
    X = (X * sign(exp)) % p

    return X, Z