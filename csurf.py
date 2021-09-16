from framework import *
exec("from src.gae_%s import *" % setting.style)

print("")
print(f'Exponent of two:\t{exponent_of_two}')
print(f'Extra cofactor 3:\t{cofactor_3}')
print(f'Extra cofactor 5:\t{cofactor_5}')
print(f'Extra cofactor 7:\t{cofactor_7}')
print("")

assert(exponent_of_two >= 3)

from src.radical.crad_4 import act_with_four_on_Montgomery, act_with_four_on_Montgomery_pro, e_2
from src.radical.crad_3 import e_3
#from src.radical.crad_5 import e_5
#from src.radical.crad_7 import e_7

''' -------------------------------------------------------------------------------------
    Number of degree-(l_i) isogeny constructions to be performed: m_i
    ------------------------------------------------------------------------------------- '''

# ==========================================================================
# Reading the suitable bounds
raw = '_radicals' * setting.radicals
vectorbound = "csurf_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + raw
exec("from tmp.%s import m, e_2, e_3, e_5, e_7" % vectorbound)
# ==========================================================================
temporal_m = list(set(m))
if len(temporal_m) > 1:
    # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
    LABEL_m = 'different_bounds'
else:
    # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
    LABEL_m = 'with_same_bounds'

if setting.verbose:
    tunned = '-tunned'
else:
    tunned = ''

try:

    # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
    m_prime = [ geometric_serie(m[k], L[k]) for k in range(n) ]
    r_out, L_out, R_out = rounds(m_prime[::-1], n)
    for j in range(0, len(r_out), 1):

        R_out[j] = list([L[::-1][k] for k in R_out[j]])
        L_out[j] = list([L[::-1][k] for k in L_out[j]])

    f = open('./strategies/' + 'csurf-' + setting.prime  + '-' + setting.style + '-' + 'hvelu' + '-' + tunned + '-' + LABEL_m + '-e' + str(setting.exponent) + raw.replace('_','-'))
    print("// Reading Strategies from the file ./strategies/csurf-%s-%s-%s-%s-%s-e%s%s" % (setting.prime, setting.style, 'hvelu', tunned, LABEL_m, str(setting.exponent), raw.replace('_','-')) )
    S_out = []
    for i in range(0, len(r_out), 1):

        tmp = f.readline()
        tmp = [ int(b) for b in tmp.split() ]
        S_out.append(tmp)

    f.close()

except IOError:

    print("// Constructing Strategies")
    C_out, L_out, R_out, S_out, r_out = strategy_block_cost(L[::-1], m[::-1])
    print("// Storing Strategies into the file ./strategies/csurf-%s-%s-%s-%s-%s-e%s%s" % (setting.prime, setting.style, 'hvelu', tunned, LABEL_m, str(setting.exponent), raw.replace('_','-')) )
    f = open('./strategies/' + 'csurf-' + setting.prime  + '-' + setting.style + '-' + 'hvelu' + '-' + tunned + '-' + LABEL_m + '-e' + str(setting.exponent) + raw.replace('_','-'),'w')
    for i in range(0, len(r_out)):

        f.writelines(' '.join([ str(tmp) for tmp in S_out[i]]) + '\n')

    f.close()

print("// All the experiments are assuming S = %1.6f x M and a = %1.6f x M. The measures are given in millions of field operations." % (SQR, ADD))

data_crads = {'S':[]}
if setting.algorithm == 'csurf' and setting.radicals:
    #m = [e_3] + [e_5] + [e_7] + m[3:]
    m = [e_3] + m[1:]
    for i in range(len(S_out)):
        # It doesn't matter the last small ell_i, it only add log2(n_i) scalar point multiplications
        TMP, _ = dynamic_programming_algorithm(L_out[i] + [3], len(L_out[i]) + 1)
        data_crads['S'].append(TMP)

def validate(pub : int):
    return issupersingular([ fp_add(pub, 2), 4])

# =================================
# Public and Private key generation
def keygen():

    # Recall, the affine Montgomery curve coefficient A = 0 has a projective representation as (A + 2 : 4) = (2 : 4)
    priv = [random.randint(-e_2, e_2)] + random_key(m)

    set_zero_ops()
    if setting.radicals:
        if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
            # This branch is focused when m degree-ell isogeny constructions are required for each ell playing on the GAE
            #pub = GAE([2, 4], priv[1:], [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m, crads = {**data_crads, 'L':[3, 5, 7]})
            pub = GAE([2, 4], priv[1:], [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m, crads = {**data_crads, 'L':[3]})
        else:
            # This branch is centered when m_i degree-(ell_i) isogeny constructions are required per each ell playing on the GAE
            #pub = GAE([2, 4], priv[1:], L_out, R_out, S_out, r_out, m, crads = {**data_crads, 'L':[3, 5, 7]})
            pub = GAE([2, 4], priv[1:], L_out, R_out, S_out, r_out, m, crads = {**data_crads, 'L':[3]})
    else:
        if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
            # This branch is focused when m degree-ell isogeny constructions are required for each ell playing on the GAE
            pub = GAE([2, 4], priv[1:], [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m)
        else:
            # This branch is centered when m_i degree-(ell_i) isogeny constructions are required per each ell playing on the GAE
            pub = GAE([2, 4], priv[1:], L_out, R_out, S_out, r_out, m)
    
    #rewrite to proj_coeffs used in rad isogenies
    X = fp_add(pub[0], pub[0])
    Z = pub[1]
    X = fp_sub(X, Z)
    X = fp_add(X, X)
    
    cost_velu = get_ops()

    set_zero_ops()
    # Degree-2 radical isogeny chain on the floor --> surface --> floor
    #pub = act_with_two_on_Montgomery(pub, priv[0])
    
    # Degree-4 radical isogeny chain on the floor --> surface --> floor
    pub = act_with_four_on_Montgomery_pro(X, Z, priv[0])

    cost_rad = get_ops()
    cost = [x + y for x, y in zip(cost_velu, cost_rad)]
    
    return priv, pub, cost_velu, cost_rad, cost

# ==========================
# Secret Sharing Computation
def derive(priv, pub : int):

    # Ensuring the input public curve is supersingular
    assert(validate(pub))

    # Recall, the affine Montgomery curve coefficient A has a projective representation as (A + 2 : 4)
    curve = [ fp_add(pub, 2), 4]

    set_zero_ops()
    if setting.radicals:
        if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
            # This branch is focused when m degree-ell isogeny constructions are required for each ell playing on the GAE
            #ss = GAE(curve, priv[1:], [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m, crads = {**data_crads, 'L':[3, 5, 7]})
            ss = GAE(curve, priv[1:], [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m, crads = {**data_crads, 'L':[3]})
        else:
            # This branch is centered when m_i degree-(ell_i) isogeny constructions are required per each ell playing on the GAE
            #ss = GAE(curve, priv[1:], L_out, R_out, S_out, r_out, m, crads = {**data_crads, 'L':[3, 5, 7]})
            ss = GAE(curve, priv[1:], L_out, R_out, S_out, r_out, m, crads = {**data_crads, 'L':[3]})
    else:
        if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
            # This branch is focused when m degree-ell isogeny constructions are required for each ell playing on the GAE
            ss = GAE(curve, priv[1:], [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m)
        else:
            # This branch is centered when m_i degree-(ell_i) isogeny constructions are required per each ell playing on the GAE
            ss = GAE(curve, priv[1:], L_out, R_out, S_out, r_out, m)        

    #rewrite to proj_coeffs used in rad isogenies
    X = fp_add(ss[0], ss[0])
    Z = ss[1]
    X = fp_sub(X, Z)
    X = fp_add(X, X)

    cost_velu = get_ops()

    set_zero_ops()
    # Degree-2 radical isogeny chain on the floor --> surface --> floor
    #ss = act_with_two_on_Montgomery(ss, priv[0])

    # Degree-4 radical isogeny chain on the floor --> surface --> floor
    ss = act_with_four_on_Montgomery_pro(X, Z, priv[0])

    cost_rad = get_ops()    
    cost = [x + y for x, y in zip(cost_velu, cost_rad)]
    
    return ss, cost_velu, cost_rad, cost

# =======================================
def pubkey_print(label : str, pub : int):
    print("%s := 0x%X;" % (label, pub))
    return None

# ========================================
def privkey_print(label : str, priv):
    print("%s := ( %s );" % (label, ' '.join(map(str, priv))) )
    return None