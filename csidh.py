from framework import *
exec("from src.gae_%s import *" % setting.style)

print("")
print(f'Exponent of two:\t{exponent_of_two}')
print(f'Extra cofactor 3:\t{cofactor_3}')
print(f'Extra cofactor 5:\t{cofactor_5}')
print(f'Extra cofactor 7:\t{cofactor_7}')
print("")

''' -------------------------------------------------------------------------------------
    Number of degree-(l_i) isogeny constructions to be performed: m_i
    ------------------------------------------------------------------------------------- '''

# ==========================================================================
# Reading the suitable bounds
vectorbound = "csidh_" + setting.prime  + "_" + setting.style + "_m" + str(str(setting.exponent))
exec("from tmp.%s import m" % vectorbound)

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

    f = open('./strategies/' + 'csidh-' + setting.prime  + '-' + setting.style + '-' + 'hvelu' + '-' + tunned + '-' + LABEL_m + '-e' + str(setting.exponent))
    print("// Reading Strategies from the file ./strategies/csidh-%s-%s-%s-%s-%s-e%s" % (setting.prime, setting.style, 'hvelu', tunned, LABEL_m, str(setting.exponent)) )
    S_out = []
    for i in range(0, len(r_out), 1):

        tmp = f.readline()
        tmp = [ int(b) for b in tmp.split() ]
        S_out.append(tmp)

    f.close()

except IOError:

    print("// Constructing Strategies")
    C_out, L_out, R_out, S_out, r_out = strategy_block_cost(L[::-1], m[::-1])
    print("// Storing Strategies into the file ./strategies/csidh-%s-%s-%s-%s-%s-e%s" % (setting.prime, setting.style, 'hvelu', tunned, LABEL_m, str(setting.exponent)) )
    f = open('./strategies/' + 'csidh-' + setting.prime  + '-' + setting.style + '-' + 'hvelu' + '-' + tunned + '-' + LABEL_m + '-e' + str(setting.exponent),'w')
    for i in range(0, len(r_out)):

        f.writelines(' '.join([ str(tmp) for tmp in S_out[i]]) + '\n')

    f.close()

print("// All the experiments are assuming S = %1.6f x M and a = %1.6f x M. The measures are given in millions of field operations." % (SQR, ADD))

def validate(pub : int):
    return issupersingular([ fp_add(pub, 2), 4])

# =================================
# Public and Private key generation
def keygen():

    # Recall, the affine Montgomery curve coefficient A = 0 has a projective representation as (A + 2 : 4) = (2 : 4)
    priv = random_key(m)

    set_zero_ops()
    if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
        # This branch is focused when m degree-ell isogeny constructions are required for each ell playing on the GAE
        pub = GAE([2, 4], priv, [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m)
    else:
        # This branch is centered when m_i degree-(ell_i) isogeny constructions are required per each ell playing on the GAE
        pub = GAE([2, 4], priv, L_out, R_out, S_out, r_out, m)

    pub = coeff(pub)
    cost = get_ops()

    cost_velu = cost
    cost_rad = [0,0,0]
    
    return priv, pub, cost_velu, cost_rad, cost

# ==========================
# Secret Sharing Computation
def derive(priv, pub : int):

    # Ensuring the input public curve is supersingular
    assert(validate(pub))

    # Recall, the affine Montgomery curve coefficient A has a projective representation as (A + 2 : 4) = (2 : 4)
    curve = [ fp_add(pub, 2), 4]

    set_zero_ops()
    if (len(temporal_m) == 1) or ((len(temporal_m) == 2) and (0 in temporal_m)):
        # This branch is focused when m degree-ell isogeny constructions are required for each ell playing on the GAE
        ss = GAE(curve, priv, [L_out[0]], [R_out[0]], [S_out[0]], [temporal_m[-1]], m)
    else:
        # This branch is centered when m_i degree-(ell_i) isogeny constructions are required per each ell playing on the GAE
        ss = GAE(curve, priv, L_out, R_out, S_out, r_out, m)

    ss = coeff(ss)
    cost = get_ops()
    
    cost_velu = cost
    cost_rad = [0,0,0]
    
    return ss, cost_velu, cost_rad, cost

# =======================================
def pubkey_print(label : str, pub : int):
    print("%s := 0x%X;" % (label, pub))
    return None

# ========================================
def privkey_print(label : str, priv):
    print("%s := ( %s );" % (label, ' '.join(map(str, priv))) )
    return None