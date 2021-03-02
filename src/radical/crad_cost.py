# -*- coding: utf-8 -*-
"""
All cost is measured in multiplications (M).
C_S is a global coefficient for the cost of squarings (S) in terms of M
So usually C_S = 0.8 or C_S = 1.0
C_A is a global coefficient for the cost of additions (A) in terms of M
So usually C_A = 0.001 or C_A = 0.01

We assume that with addition chains, an exponentiation can be close to (1+C_exp)*log(N) mults
where C_exp is 0.2 (instead of the usual 0.5 with square and multiply)
In practice, C_exp is usually as low as 0.1 or 0.05 for large N

For specific primes, addchains can be used to calculate the exact number of field operations
required for the required exponents. 
For the primes used in the paper, these short addition chains can be found in cradprimes folder
"""
from math import log, ceil
from src.fp import p
from src.montgomery import C_xMUL, global_L, measure
from src.radical.crad_general import sq_exp, tri_exp, quart_exp, quint_exp, sept_exp, novem_exp

C_S = 1.0
C_A = 0.01
C_exp = 0.2

def cost_exponentiation(N):
    sqrs = log(N, 2)
    mult = C_exp*log(N,2)
    
    cost_addchain = C_S*sqrs + mult
    
    return cost_addchain

cost_inv = cost_exponentiation(p-2)
cost_sq_root = cost_exponentiation(sq_exp)
cost_tri_root = cost_exponentiation(tri_exp)
cost_quart_root = cost_exponentiation(quart_exp)
cost_quint_root = cost_exponentiation(quint_exp)
cost_sept_root = cost_exponentiation(sept_exp)
cost_novem_root = cost_exponentiation(novem_exp)


def cost_rad(ell, exp):
    isogeny_cost = exp*cost_single_isogeny(ell)
    overhead = cost_overhead(ell)
    if ell % 2 == 0:
        sampling = 0
    else:
        sampling = cost_sampling(ell)
    
    total = isogeny_cost + overhead + sampling
    
    return ceil(total)

#all of these costs have been calculated from the other crad_#.py files
def cost_single_isogeny(ell):
    
    if ell == 2:
        cost = cost_sq_root + C_S*1 + 1 + C_A*4
    elif ell == 3:
        cost = cost_tri_root + C_S*0 + 2 + C_A*10
    elif ell == 4:
        cost = cost_quart_root + C_S*5 + 3 + C_A*3
    elif ell == 5:
        cost = cost_quint_root + C_S*6 + 8 + C_A*8
    elif ell == 7:
        cost = cost_sept_root + C_S*3 + 28 + C_A*20
    elif ell == 9:
        cost = cost_novem_root + cost_inv + C_S*5 + 35 + C_A*91        

    return cost

cost_mont_min = 1*cost_sq_root + 0*cost_tri_root + 1*cost_inv + C_S*1 + 1 + C_A*2
cost_min_mont = 1*cost_sq_root + 0*cost_tri_root + 1*cost_inv + C_S*1 + 1 + C_A*2
cost_mont_min_mont = cost_mont_min + cost_min_mont

cost_min_tate4 = 0*cost_sq_root + 0*cost_tri_root + 1*cost_inv + C_S*1 + 7 + C_A*10
cost_tate4_min = 2*cost_sq_root + 0*cost_tri_root + 2*cost_inv + C_S*1 + 3 + C_A*7
cost_min_tate4_min = cost_min_tate4 + cost_tate4_min

cost_mont_tate = 1*cost_sq_root + 0*cost_tri_root + 1*cost_inv + C_S*1 + 9 + C_A*10
cost_tate_mont = 2*cost_sq_root + 1*cost_tri_root + 1*cost_inv + C_S*7 + 23 + C_A*49
cost_mont_tate_mont = cost_mont_tate + cost_tate_mont

#all of these costs have been calculated from the other crad_#.py files
def cost_overhead(ell):
    
    if ell == 2:
        cost_ops = 3*cost_sq_root + 2*cost_inv + C_S*2 + 6 + C_A*12
    elif ell == 3:
        cost_ops = 0*cost_sq_root + 0*cost_inv + C_S*0 + 2 + C_A*0
    elif ell == 4:
        cost_ops = 3*cost_sq_root + 1*cost_inv + C_S*5 + 6 + C_A*7
    elif ell == 5:
        cost_ops = 0*cost_sq_root + 1*cost_inv + C_S*4 + 5 + C_A*1
    elif ell == 7:
        cost_ops = 0*cost_sq_root + 1*cost_inv + C_S*1 + 3 + C_A*4
    elif ell == 9:
        cost_ops = 0*cost_sq_root + 1*cost_inv + C_S*3 + 4 + C_A*9
        
    if ell == 2:
        cost_functs = cost_mont_min_mont
    elif ell == 3:
        cost_functs = cost_mont_tate_mont
    elif ell == 4:
        cost_functs = cost_mont_min_mont + cost_min_tate4_min
    elif ell == 5:
        cost_functs = cost_mont_tate_mont
    elif ell == 7:
        cost_functs = cost_mont_tate_mont
    elif ell == 9:
        cost_functs = cost_mont_tate_mont
    
    cost = cost_ops + cost_functs
    return cost


def cost_sampling(ell):

    if ell != 2 and ell != 4:
        # Let's assume the best scenario, only one random sampling of points is performed.
        # Thus, the cost of a  sampling is just the cost of [(p+1) / ell]P
        return measure(sum(C_xMUL) - C_xMUL[global_L.index(ell)])
    else:
        # Radical isogenies of degree 2 and 4 do not require sampling!
        return 0