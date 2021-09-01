from math import copysign
from src.fp import *
from src.montgomery import xDBL, xMUL

raw = '_radicals' * setting.radicals
vectorbound = "csurf_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + raw
exec("from tmp.%s import m, e_2, e_3, e_5, e_7" % vectorbound)

# ---

sq_exp = (p + 1)//4             #equiv to taking square root
quart_exp = (p + 1)//8          #if p suitable, equiv to taking fourth root
tri_exp = (2*p - 1)//3          #equiv to taking third root
quint_exp = (3*p-2)//5          #equiv to taking fifth root
sept_exp = (4*p-3)//7           #equiv to taking seventh root
novem_exp = (5*p - 4)//9        #equiv to taking ninth root

# --- Radical exponents: procedure named as crad_###()
rnd = random.randint(2, p)
exec("from radical.inv_%s import inv_%s as crad_inv" % (setting.prime[1:], setting.prime[1:]))                          # inv_###
label = 'inv'
assert(fp_inv(rnd) == crad_inv(rnd))
#print(f"// Tested fp_{label}(###) vs. crad_{label}(###)")
for label in ['sq', 'tri', 'quart', 'quint', 'sept', 'novem']:
    exec("from radical.%s_%s import %s_%s as crad_%s" % (label, setting.prime[1:], label, setting.prime[1:], label))    # label_###
    exec("fp_exp(%s, %s_exp) == crad_%s(%s)" % (str(rnd),label, label,str(rnd)))
    #print(f"// Tested fp_exp(###, {label}_exp) vs. crad_{label}(###)")

#print("\nEverything is okay!\n")

inv_2 = crad_inv(2)       # 1 / 2
inv_3 = crad_inv(3)       # 1 / 3
inv_4 = fp_sqr(inv_2)   # 1 / 4


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
    inv = crad_sq(inv)            #(A^2 + 4)^sq_exp
    inv = crad_inv(inv)             #1/(A^2 + 4)^sq_exp
    output = fp_mul(output,inv)   #(-2*A)/((A^2 + 4)^sq_exp)
    return output


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
    inv = crad_sq(inv)      

    return output, inv


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
    inv = crad_sq(inv)     #(4 - A^2)^sq_exp
    inv = crad_inv(inv)             #1/((4 - A^2)^sq_exp)
    output = fp_mul(output,inv)   #(-2*A)/((4 - A^2)^sq_exp)
    return output


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
    inv = crad_sq(inv)      
    
    return output, inv


def Montgomery_to_Tate(A, xP):
    '''
    From the Radical isogenies code: 
    This functions transforms
        FROM: a Montgomery coefficient A on the floor
                    y^2 = x^3 + A*x^2 + x
        TO: the coefficients of a Tate normal form of X1(4)
                y^2 + M*xy + N*y = x^3 + N*x^2
    obtained by an isomorphism that puts a point with x-coordinate xP at (0,0), 
    where xP is assumed to be an ell-torsion point's x-coordinate. 
    The choice of sign for the y-coordinate corresponding to xP is arbitrary but fixed. 
    Note that in classical notation, M = 1-c and N = -b.
    The specific case for degree 3 can be found in crad_3.py
    We use the general naming conventions [u,r,s,t] for Weierstrass-isomorphisms.
'''

    r = xP
    tmp = fp_add(A,r)                   #A+r
    t = fp_mul(r,tmp)                   #r*(A+r)
    t = fp_add(1, t)                    #(1 + r*(A+r))
    t = fp_mul(r,t)                     #r * (1 + r*(A+r))
    t = crad_sq(t)                #(r*(1+r*(A+r)))^sq_exp
    twot_inv = fp_add(t,t)              #2*t
    twot_inv = crad_inv(twot_inv)         #1/(2*t)
    s = fp_add(tmp, tmp)                #2*A + 2*r
    s = fp_add(s, r)                    #2*A + 3*r
    s = fp_mul(r,s)                     #r*(2*A + 3*r)
    s = fp_add(1,s)                     #1 + r*(2*A + 3*r)
    s = fp_mul(s, twot_inv)             #(1 + r*(2*A + 3*r))/(2*t)
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
    twot_inv = crad_inv(twot_inv)         #1/(2*t)
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
    B = crad_inv(B)                           #1/4*A
    B = fp_sub(B,2)                         #1/4*A - 2
    eps = fp_sqr(B)                         #B^2
    eps = fp_sub(eps, 4)                    #B^2 - 4
    eps = crad_sq(eps)                      #(B^2 - 4)^sq_exp
    output = fp_mul(3,eps)                  #3*eps
    output = fp_add(-B, output)             #(-B + 3*eps)
    out_inv = fp_sub(B,eps)                 #(B-eps)
    out_inv = fp_mul(eps, out_inv)          #eps*(B-eps)
    out_inv = fp_add(out_inv,out_inv)       #eps*(B-eps)*2
    out_inv = crad_sq(out_inv)              #(eps*(B-eps)*2)^sq_exp
    out_inv = crad_inv(out_inv)             #1/(#(eps*(B-eps)*2)^sq_exp)
    output = fp_mul(output, out_inv)        #(-B + 3*eps)/(#(eps*(B-eps)*2)^sq_exp)
    
    return output

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
    eps = crad_sq(eps)               #((Bp + 8*X)(Bp - 8*X)*(4X)^2)^sq_exp

    _3eps = fp_add(eps, fp_add(eps, eps))
    ttmp = fp_mul(Bp, _4AX)
    output = fp_sub(_3eps, ttmp)            #-Bp*(4x) + 3*eps

    inv = fp_mul(Bp, _4AX)
    inv = fp_sub(inv, eps)
    inv = fp_mul(eps, inv)
    inv = fp_add(inv, inv)
    inv = crad_sq(inv)               #(2*eps*(BP(4X) - eps))^sq_exp
    
    return output, inv

def Weier_to_Montgomery(coeffs):
    '''
   This function takes a 5-tuple [a1,a2,a3,a4,a6] of a (long) Weierstrass form and
   calculates the corresponding Montgomery coefficient.
   Note that not all curves have a Montgomery representation, and we just assume
   the input does have one.
   Furthermore, even if the curve has a Montgomery form, the formulas only work for
   the case where the rational 2-torsion of the curve is Z/2Z. The reason is that
   this allows us to solve a cubic polynomial without using field extensions (i.e.
   we need to avoid casus irreducibilis for the cubic polynomial). In particular,
   for CSIDH, this means the computations need to be done on curves on the floor.
   '''
    a1 = coeffs[0]
    a2 = coeffs[1]
    a3 = coeffs[2]
    a4 = coeffs[3]
    a6 = coeffs[4]

    s = fp_sub(0,a1)
    s = fp_mul(s, inv_2)             #-a1/2

    B = fp_sqr(a1)
    B = fp_mul(B, inv_4)
    B = fp_add(B, a2)                #(a1^2/4 + a2)

    C = fp_mul(a1, a3)
    C = fp_mul(C, inv_2)
    C = fp_add(C, a4)                #(a1*a3/2 + a4)

    D = fp_sqr(a3)
    D = fp_mul(D, inv_4)
    D = fp_add(D, a6)                #a3^2/4 + a6

    d_tmp = fp_sqr(B)
    d0 = fp_sub(d_tmp, fp_mul(3, C))   #B^2-3*C

    d1 = fp_mul(2, d_tmp)
    d1 = fp_sub(d1, fp_mul(9, C))
    d1 = fp_mul(B, d1)
    d1 = fp_add(d1, fp_mul(27,D))       #B*(2*B^2-9*C)+27*D   

    d0c = fp_sqr(d0)
    d0c = fp_mul(d0c, d0)
    sq_root = fp_sqr(d1)
    sq_root = fp_sub(sq_root, fp_mul(4, d0c))
    sq_root = crad_sq(sq_root)          #(d1^2-4*d0^3)^sq_exp

    cb_root = fp_add(d1, sq_root)
    cb_root = fp_mul(cb_root, inv_2)
    cb_root = crad_tri(cb_root)         #((d1+sq_root)/2)^tri_exp;

    r = crad_inv(cb_root)
    r = fp_mul(d0, r)
    r = fp_add(cb_root, r)
    r = fp_add(B, r)
    r = fp_mul(-r, inv_3)               #-(B+cuberoot+delta0/cuberoot)/3

    t = fp_mul(r, a1)
    t = fp_add(a3, t)
    t = fp_mul(-t, inv_2)               #-(a3+r*a1)/2

    u2_2 = fp_mul(s, a3)
    u2_3 = fp_mul(r, a2)
    u2_3 = fp_mul(2, u2_3)
    u2_4 = fp_mul(r,s)
    u2_4 = fp_add(t, u2_4)
    u2_4 = fp_mul(u2_4, a1)
    u2_5 = fp_sqr(r)
    u2_5 = fp_mul(3, u2_5)
    u2_6 = fp_mul(s,t)
    u2_6 = fp_add(u2_6, u2_6)

    u2 = fp_sub(a4, u2_2)
    u2 = fp_add(u2, u2_3)
    u2 = fp_sub(u2, u2_4)
    u2 = fp_add(u2, u2_5)
    u2 = fp_sub(u2, u2_6)
    u2 = crad_sq(u2)                    #(a4 -s*a3 +2*r*a2 -(t+r*s)*a1 +3*r^2 -2*s*t)^sq_exp

    output_2 = fp_mul(s,a1)
    output_3 = fp_mul(3,r)
    output_4 = fp_sqr(s)
    output = fp_sub(a2, output_2)
    output = fp_add(output, output_3)
    output = fp_sub(output, output_4)
    output = fp_mul(output, crad_inv(u2))         #(a2-s*a1+3*r-s^2)/u2

    return output

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


# Next function correspond with the affine version of elligator
def sampling_ell_order_point(A, ell):

    u = random.randint(2, p_minus_one_halves)
    u_squared = fp_sqr(u)

    u_squared_plus_one  = fp_add(u_squared, 1)
    u_squared_minus_one = fp_sub(u_squared, 1)

    A_times_u_squared_minus_one = fp_mul(A, u_squared_minus_one)

    tmp = fp_sqr(A)
    tmp = fp_mul(tmp, u_squared)
    aux = fp_sqr(u_squared_minus_one)
    tmp = fp_add(tmp, aux)
    tmp = fp_mul(A_times_u_squared_minus_one, tmp)

    alpha, beta = 0, u
    alpha, beta = fp_cswap(alpha, beta, tmp == 0)
    u_squared_plus_one = fp_mul(alpha, u_squared_plus_one)
    alpha = fp_mul(alpha, u_squared_minus_one)

    Tp_X = fp_add(A, alpha)
    Tm_X = fp_mul(A, u_squared)
    Tm_X = fp_add(Tm_X, alpha)
    Tm_X = fp_sub(0, Tm_X)

    tmp = fp_add(tmp, u_squared_plus_one)
    Tp_X, Tm_X = fp_cswap(Tp_X, Tm_X, (1 - jacobi(tmp, p)) // 2 )

    Tp_proj = [Tp_X, u_squared_minus_one]
    Tm_proj = [Tm_X, u_squared_minus_one]

    # Next, make Tp and Tm of torsion (p+1)/2^e
    A_proj = [fp_add(A, 2), 4]
    for i in range(0, exponent_of_two, 1):
        Tp_proj = xDBL(Tp_proj, A_proj)
        Tm_proj = xDBL(Tm_proj, A_proj)

    # Finally, we make Tp and Tm of torsion (p+1)/(2^e * the product of L[ell]' != L[ell])
    for i in range(0, n, 1):
        if i != ell:
            # Multipliying by L[i] with L[i] != L[ell]
            Tp_proj = xMUL(Tp_proj, A_proj, i)
            Tm_proj = xMUL(Tm_proj, A_proj, i)

    # Making torsion-L[ell] for L[ell] != 3.
    if ell != 0 and cofactor_3:
        # Eliminating one extra factor 3
        assert(L[0] == 3)
        Tp_proj = xMUL(Tp_proj, A_proj, 0)
        Tm_proj = xMUL(Tm_proj, A_proj, 0)

    # Making torsion-L[ell] for L[ell] != 5.
    if ell != 1 and cofactor_5:
        # Eliminating one extra factor 5
        assert(L[1] == 5)
        Tp_proj = xMUL(Tp_proj, A_proj, 1)
        Tm_proj = xMUL(Tm_proj, A_proj, 1)

    # Making torsion-L[ell] for L[ell] != 7.
    if ell != 2 and cofactor_7:
        # Eliminating one extra factor 7
        assert(L[2] == 7)
        Tp_proj = xMUL(Tp_proj, A_proj, 2)
        Tm_proj = xMUL(Tm_proj, A_proj, 2)

    # Otherwise, we return either 
    # - a pair of torsion-9 points. Recall torsion-9 != 9-order, e.g., a 3-order point is a torsion-9 point
    # - a pair of torsion-25 points. Recall torsion-25 != 25-order, e.g., a 5-order point is a torsion-25 point
    # - a pair of torsion-49 points. Recall torsion-49 != 49-order, e.g., a 7-order point is a torsion-49 point

    return Tp_proj, Tm_proj