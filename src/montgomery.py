# Loading the library corresponding to the arithmetic in F_p
from src.fp import *

bitlength = lambda x: len(bin(x)[2:])                       # number of bits
hamming_weight = lambda x: bin(x).count("1");               # hamming weight: number of bits equal 1

'''
    dacs()
    inputs: a small odd prime number l, three integer numbers, and a list
    output: all the differential additions chains corresponding with the input l

    NOTE: this is a recursive approach
'''
def dacs(l, r0, r1, r2, chain):

    if r2 == l:

        return [(chain, r2)]
    elif r2 < l and len(chain) <= 1.5*log(l,2):

        return dacs(l, r0, r2, r2 + r0, chain + [1]) + dacs(l, r1, r2, r2 + r1, chain + [0])
    else:
        return []


'''
    sdac()
    input: a small odd prime number l
    output: the shortest differential additions chains corresponding with the input l

    NOTE: this function uses a recursive function
'''
def sdac(l):

    all_dacs = dacs(l, 1, 2, 3, [])
    return min(all_dacs, key=lambda t: len(t[0]))[0]

# -------------------------------------------------------------------------------------------------------------------------------
#########################################################################################################
SQR = 1.00;                                                 # In F_p, we have SQR_{F_p} = SQR x MUL_{F_p}
ADD = 0.00;                                                 # In F_p, we have ADD_{F_p} = ADD x MUL_{F_p}
#########################################################################################################

measure = lambda x: (x[0] + SQR * x[1] + ADD * x[2])

try:

    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
    f = open('./sdacs/' + setting.prime)
    if( (sys.argv[0] != 'header.py')  and (sys.argv[0] != 'tunned-parameters.py')):
        print("// Reading Shortest Differential Addition Chains (SDACs) from the file ./sdacs/%s" % setting.prime)

    SDACS = []
    SDACS_LENGTH = []
    for i in range(0, n, 1):

        tmp = f.readline()
        tmp = [ int(b) for b in tmp.split() ]
        SDACS.append(tmp)
        SDACS_LENGTH.append(len(tmp))

    f.close()

except IOError:

    if( (sys.argv[0] != 'header.py')  and (sys.argv[0] != 'tunned-parameters.py')):
        print("// Computing Shortest Differential Addition Chains (SDACs)")

    SDACS = list(map(sdac, L))                              # Shortest Differential Addition Chains for each small odd prime l in L
    SDACS_LENGTH = [ len(sdacs) for sdacs in SDACS]

    # Storing the SDACS to be used
    if( (sys.argv[0] != 'header.py')  and (sys.argv[0] != 'tunned-parameters.py')):
        print("// Storing Shortest Differential Addition Chains (SDACs) into the file ./sdacs/%s" % setting.prime)

    f = open('./sdacs/' + setting.prime,'w')
    for i in range(0, n, 1):

        f.writelines(' '.join([ str(tmp) for tmp in SDACS[i]]) + '\n')
    f.close()

global_L = list(L)
cMUL  = lambda l: numpy.array( [ 4.0 * (SDACS_LENGTH[L.index(l)] + 2), 2.0 * (SDACS_LENGTH[L.index(l)] + 2), 6.0 * (SDACS_LENGTH[L.index(l)] + 2) - 2.0] )
C_xMUL  = list(map(cMUL,  global_L))   # list of the costs of each [l]P

''' ----------------------------------------------------------------------
    coeff()
    input : projective Montgomery constants A24 := A + 2C and C24 := 4C
            where E : y^2 = x^3 + (A/C)*x^2 + x
    output: the affine Montgomery coefficient A/C
    ---------------------------------------------------------------------- '''
def coeff(A):
    output = fp_add(A[0], A[0])         # (2 * A24)
    output = fp_sub(output, A[1])       # (2 * A24) - C24
    C24_inv = fp_inv(A[1])              # 1 / (C24)
    output = fp_add(output, output)     # 4*A = 2[(2 * A24) - C24]
    output = fp_mul(output, C24_inv)    # A/C = 2[(2 * A24) - C24] / C24

    return output

# isinfinity(P) determines if x(P) := (XP : ZP) = (1 : 0)
def isinfinity(P):
    return P[1] == 0

# areequal(P, Q) determines if x(P) = x(Q)
def areequal(P, Q):
    return fp_mul(P[0], Q[1]) == fp_mul(P[1], Q[0])

''' ----------------------------------------------------------------------
    xDBL()
    input : a projective Montgomery x-coordinate point x(P) := XP/ZP, and
            the  projective Montgomery constants A24:= A + 2C and C24:=4C 
            where E : y^2 = x^3 + (A/C)*x^2 + x
    output: the projective Montgomery x-coordinate point x([2]P)
    ---------------------------------------------------------------------- '''
def xDBL(P, A):

    t_0 = fp_sub(P[0],P[1])
    t_1 = fp_add(P[0],P[1])
    t_0 = fp_sqr(t_0)
    t_1 = fp_sqr(t_1)
    Z = fp_mul(A[1], t_0);
    X = fp_mul(Z, t_1);
    t_1 = fp_sub(t_1, t_0);
    t_0 = fp_mul(A[0], t_1);
    Z = fp_add(Z, t_0);
    Z = fp_mul(Z, t_1);

    return [X, Z]

''' ----------------------------------------------------------------------
    xADD()
    input : the projective Montgomery x-coordinate points x(P) := XP/ZP, 
            x(Q) := XQ/ZQ, and x(P-Q) := XPQ/ZPQ
    output: the projective Montgomery x-coordinate point x(P+Q)
    ---------------------------------------------------------------------- '''
def xADD(P, Q, PQ):

    a = fp_add(P[0], P[1])
    b = fp_sub(P[0], P[1])
    c = fp_add(Q[0], Q[1])
    d = fp_sub(Q[0], Q[1])
    a = fp_mul(a, d)
    b = fp_mul(b, c)
    c = fp_add(a, b)
    d = fp_sub(a, b)
    c = fp_sqr(c)
    d = fp_sqr(d)
    X = fp_mul(PQ[1], c)
    Z = fp_mul(PQ[0], d)
    return [X, Z]

''' ----------------------------------------------------------------------
    xMUL()
    input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
            projective Montgomery constants A24:= A + 2C and C24:=4C where 
            E : y^2 = x^3 + (A/C)*x^2 + x, and an positive integer j
    output: the projective Montgomery x-coordinate point x([L[j]]P)
    ---------------------------------------------------------------------- '''
    # Modificar esta parte para usar cadenas de addicion
def xMUL(P, A, j):

    P2= xDBL(P, A)
    R = [P, P2, xADD(P2, P, P)]

    for i in range(SDACS_LENGTH[j] - 1, -1, -1):

        if isinfinity(R[SDACS[j][i]]):
            T = xDBL(R[2], A)
        else:
            T = xADD(R[2], R[SDACS[j][i] ^ 1], R[SDACS[j][i]])

        R[0] = list(R[SDACS[j][i] ^ 1])
        R[1] = list(R[2])
        R[2] = list(T)

    return R[2]

''' ----------------------------------------------------------------------
    prime_factors()
    input : a projective Montgomery x-coordinate point x(P) := XP/ZP, the
            projective Montgomery constants A24:= A + 2C and C24:=4C where 
            E : y^2 = x^3 + (A/C)*x^2 + x, and subset of |[0, n]|
    output: the projective Montgomery x-coordinate points x([(p+1) / l_0]P),
            x([(p+1) / l_1]P), ..., x([(p+1) / l_{n-1}]P).
    ---------------------------------------------------------------------- '''
def prime_factors(P, A, points):
    n = len(points)
    if n == 1:
        # In this recursion level we have an order-l point
        return [P]
    elif n > 0:
        # We proceed by applying a divide-and-conquer procedure
        h = n // 2
        if h > 0:

            # 1st half
            first_half = []
            second_P = P
            for j in range(h):
                second_P = xMUL(second_P, A, points[j])
                first_half.append(points[j])

            # 2nd half
            second_half = []
            first_P = P
            for j in range(h, n):
                first_P = xMUL(first_P, A, points[j])
                second_half.append(points[j])

            return prime_factors(first_P, A, first_half) + prime_factors(second_P, A, second_half)

        return []


def elligator(A):

    Ap = fp_add(A[0], A[0])
    Ap = fp_sub(Ap, A[1])
    Ap = fp_add(Ap, Ap)
    Cp = A[1];

    u = random.randint(2, p_minus_one_halves)
    u_squared = fp_sqr(u)

    u_squared_plus_one  = fp_add(u_squared, 1)
    u_squared_minus_one = fp_sub(u_squared, 1)

    C_times_u_squared_minus_one = fp_mul(Cp, u_squared_minus_one)
    AC_times_u_squared_minus_one= fp_mul(Ap, C_times_u_squared_minus_one)

    tmp = fp_sqr(Ap)
    tmp = fp_mul(tmp, u_squared)
    aux = fp_sqr(C_times_u_squared_minus_one)
    tmp = fp_add(tmp, aux)
    tmp = fp_mul(AC_times_u_squared_minus_one, tmp)

    alpha, beta = 0, u
    alpha, beta = fp_cswap(alpha, beta, tmp == 0)
    u_squared_plus_one = fp_mul(alpha, u_squared_plus_one)
    alpha = fp_mul(alpha, C_times_u_squared_minus_one)

    Tp_X = fp_add(Ap, alpha)
    Tm_X = fp_mul(Ap, u_squared)
    Tm_X = fp_add(Tm_X, alpha)
    Tm_X = fp_sub(0, Tm_X)

    tmp = fp_add(tmp, u_squared_plus_one)
    Tp_X, Tm_X = fp_cswap(Tp_X, Tm_X, (1 - jacobi(tmp, p)) // 2 )

    Tp_proj = [Tp_X, C_times_u_squared_minus_one]
    Tm_proj = [Tm_X, C_times_u_squared_minus_one]

    if cofactor_3:
        # Multiplying by 3
        assert(L[0] == 3)
        if setting.algorithm == 'csidh' or not setting.radicals:
            Tp_proj = xMUL(Tp_proj, A, 0)
            Tm_proj = xMUL(Tm_proj, A, 0)

    if cofactor_5:
        # Multiplying by 5
        assert(L[1] == 5)
        Tp_proj = xMUL(Tp_proj, A, 1)
        Tm_proj = xMUL(Tm_proj, A, 1)

    if cofactor_7:
        # Multiplying by 7
        assert(L[2] == 7)
        Tp_proj = xMUL(Tp_proj, A, 2)
        Tm_proj = xMUL(Tm_proj, A, 2)

    return Tp_proj, Tm_proj

def isfull_order(seq):
    tmp = [ not isinfinity(seq_i) for seq_i in seq ]
    return reduce(lambda x,y : (x and y), tmp)

def full_torsion_points(A):

    if setting.style != 'wd1':
        output = [ [0,0], [0,0] ]
    else:
        output = [ [0,0], [1,1] ]
    
    while [0,0] in output:

        T_p, T_m = elligator(A)
        for i in range(0, exponent_of_two, 1):
            T_p = xDBL(T_p, A)

        if isfull_order(prime_factors(T_p, A, range(0, n , 1))) and output[0] == [0,0]:
            output[0] = list(T_p)
            
        if setting.style != 'wd1':
            for i in range(0, exponent_of_two, 1):
                T_m = xDBL(T_m, A)
            if isfull_order(prime_factors(T_m, A, range(0, n , 1))) and output[1] == [0,0]:
                output[1] = list(T_m)

    return output[0], output[1]

def CrissCross(alpha, beta, gamma, delta):

    t_1 = fp_mul(alpha, delta)
    t_2 = fp_mul(beta, gamma)
    return fp_add(t_1, t_2), fp_sub(t_1, t_2)
 

def issupersingular(A):

    while(True):

        T_p, _ = elligator(A)
        if setting.algorithm == 'csurf' and setting.radicals:
            T_p = xMUL(T_p, A, 0)

        for i in range(0, exponent_of_two, 1):
            T_p = xDBL(T_p, A)
    # Removing extra factor-3
        P = prime_factors(T_p, A, range(0, n , 1))

        bits_of_the_order = 0
        for i in range(0, n, 1):
        
            if isinfinity(P[i]) == False:

                Q = xMUL(P[i], A, i)
                
                if isinfinity(Q) == False:
                    return False

                bits_of_the_order += bitlength(L[i])
                if bits_of_the_order > validation_stop:
                    return True
