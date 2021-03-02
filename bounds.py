from framework import *
exec("from src.gae_%s import *" % setting.style)
from src.radical.crad_cost import cost_rad

# -----------------------------------------------------------------
print_intv = lambda v, n: ', '.join(list(map(format, v, ['2d']*n)))

# MITM procedure
KEYSPACE = 256.00               # For ensuring 128-bits of classical security
#KEYSPACE = 384.00              # For ensuring 192-bits of classical security
# vOW Golden Collision Search
#KEYSPACE = 220.2295746338436   # For ensuring 128-bits of classical security
#KEYSPACE = 305.5629079671769   # For ensuring 192-bits of classical security

e_2 = 0
e_3 = 0
e_5 = 0
e_7 = 0
# Next function looks for the local optimal solution in a neighborhood centered at seq_i: greedy algorithm
def neighboring_intvec(seq_i, L, IN, OUT, keyspace):
    if OUT[2] >= keyspace:
        return OUT
    
    else:
        minimum = IN
        if measure(IN[1]) >= measure(OUT[1]):

            for j in seq_i:

                if OUT[0][j] > 0:            
                    current_cost, _, _, _, _ = strategy_block_cost(L, OUT[0] + basis[j])
                    tmp = neighboring_intvec(   seq_i, 
                                                L, 
                                                IN, 
                                                (OUT[0] + basis[j], current_cost, security(OUT[0] + basis[j], len(L))),
                                                keyspace
                                            )
                    if measure(minimum[1]) >= measure(tmp[1]):
                        minimum = tmp
        
        return minimum

# Looking for the suitable/optimal bounds using a greedy algorithm
def optimal_bounds(L, b, r, keyspace):

    global e_2
    global e_3
    global e_5
    global e_7

    assert(r >= 1)
    n = len(L)

    RNC, _, _, _, _ = strategy_block_cost(L, b)
    SEC = security(b, n)
    e = b

    for i in range(0, n, 1):

        # The algorithm proceed by looking the best bounds when e_i <- e_i - 1
        seq_i = [ k for k in range(n) if k != i ]
        if e[i] <= r:
            continue
        
        # Set the new possible optimal bounds
        temporal_cost, _, _, _, _ = strategy_block_cost(L, e - r*basis[i])
        (e_tmp, RNC_tmp, SEC_tmp) = neighboring_intvec( seq_i,
                                                        L, 
                                                        (e, RNC, SEC),
                                                        (e - r*basis[i], temporal_cost, security(e - r*basis[i], n)),
                                                        keyspace
                                    )

        if measure(RNC_tmp) == measure(RNC):
            # If decreasing ell_i didn't give a better running time,
            # then it is expected to be the same case (or almost neglible) for decreasing another smaller ell_j than ell_i
            # This branch is to kill local optimal convergence
            return (e, RNC)

        assert(measure(RNC_tmp) < measure(RNC))
        (e, RNC, SEC) = (e_tmp, RNC_tmp, SEC_tmp)

        print("[Security := %f]" % SEC, end="\t")
        print("decreasing: e_{" + print_intv([i], 1) + "}" +\
                ", and increasing each e_j with j != " + print_intv([i], 1) + "; current optimal running-time: %7.3f" % measure(RNC))
        print("[" + print_intv(e, n) + "]\n")
    
    # --------------------------------------------------------------------------------------------------
    return (e, RNC)

# CSURF without radical isogenies will be saved with the suffix _raw
raw = '_raw' * setting.raw

try:
    # Reading the suitable bounds
    if not setting.raw and setting.algorithm == 'csurf':
        # The use of radical isogenies assumes we have precomputed an optimal bound for CSURF using only degree-2 isogenies
        vectorbound = "csurf_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + '_raw'
    else:
        # CSURF without degree-2 isogenies takes as initial optiimal bounds the given ones by CSIDH
        vectorbound = "csidh_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent)

    exec("from tmp.%s import *" % vectorbound)
except:
    # Number of degree-(l_i) isogeny constructions to be performed: m_i
    m = setting.exponent    # Each exponent for the GAE is set as the same value "e"
    vectorbound = None

# ==========================================================================
def main():

    k = 3                   # Just a pretty way of printing the logs
    global e_2
    global e_3
    global e_5
    global e_7

    keyspace = KEYSPACE    # ---
    # Next integer vector bount is given in Onuki et al. manuscript
    print("\n_______________________________________________________________________________________________________________________________")
    print("List of small odd primes")
    printl("L", L[::-1], n // k)
    print("\nInitial integer vector of bounts (e_0, ..., e_%d)" % n)

    if setting.exponent == 1:
        if setting.style == 'wd1' or setting.style == 'df':
            assert(n >= 221)
            e = [1] * 221 + [0] * (n - 221)
        else:                 
            assert(n >= 139)
            e = [1] * 139 + [0] * (n - 139)

        # --------------------------------------------------------------------------------------------------
        f = open("./tmp/" + setting.algorithm + "_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + raw + ".py", "w")
        f.write( 'm = [' + ', '.join([ str(ei) for ei in e ]) + ']')
        if setting.algorithm == 'csurf':
            f.write( '\ne_2 = %d' % e_2)
            f.write( '\ne_3 = %d' % e_3)
            f.write( '\ne_5 = %d' % e_5)
            f.write( '\ne_7 = %d' % e_7)
        f.close()
        # --------------------------------------------------------------------------------------------------

    else:
        if vectorbound != None:
            assert(not isinstance(m, int))
            e = numpy.array(list(m)[::-1])
            printl("e", e, n // k)
            RNC, _, _, _, _ = strategy_block_cost(L[::-1], e)

            print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RNC[0] / (10.0**6), RNC[1] / (10.0**6), RNC[2] / (10.0**6), measure(RNC) / (10.0**6)) )
            print("\tSecurity ~ %f\n" % security(e, n))
            r = 1
        else:
            e = [m] * (n - 1) + [3 * m // 2]
            e = numpy.array(e)
            stop = False
            for i in range(0, n, 1):
                for j in range(0, m, 1):
                    if e[i] >= 1:
                        e = e - basis[i]
                        if security(e, n) < keyspace:
                            e = e + basis[i]
                            stop = True
                            break

                if stop:
                    break

            printl("e", e, n // k)
            RNC, _, _, _, _ = strategy_block_cost(L[::-1], e)

            print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RNC[0] / (10.0**6), RNC[1] / (10.0**6), RNC[2] / (10.0**6), measure(RNC) / (10.0**6)) )
            print("\tSecurity ~ %f\n" % security(e, n))

            print("_______________________________________________________________________________________________________________________________")
            print("We proceed by searching a better integer vector of bounds\n")
            r = 1
            for k in range(1, int(ceil( (1.0*m) / (1.0*r) ))):
                e, RNC_tmp = optimal_bounds(L[::-1], e, r, keyspace)
                if measure(RNC_tmp) == measure(RNC):
                    break
                else:
                    RNC = RNC_tmp
                # --------------------------------------------------------------------------------------------------
                f = open("./tmp/" + setting.algorithm + "_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + raw + ".py", "w")
                f.write( 'm = [' + ', '.join([ str(ei) for ei in e[::-1] ]) + ']')
                if setting.algorithm == 'csurf':
                    f.write( '\ne_2 = %d' % e_2)
                    f.write( '\ne_3 = %d' % e_3)
                    f.write( '\ne_5 = %d' % e_5)
                    f.write( '\ne_7 = %d' % e_7)
                f.close()
            print("_______________________________________________________________________________________________________________________________\n")

        if setting.algorithm == 'csurf':
            # The idea here is to apply the CSIDH-strategy but with different keyspace determined by the number of radical isogenies to be used
            e_init = e
            RNC_prev = RNC
            # In terms of cost: deg-2 < deg-3 < deg-5 < deg-7 ... Thus, we can assume e_2 >= e_3 >= e_5 >= e_7
            if not setting.raw:
                print(f'Cost assuming only degree-2 isogenies:\t{measure(RNC) + cost_rad(2, e_2)}')
                print(f'Number of degree-2 isogeny constructions on the surface:\t{e_2}')
                RNC = (measure(RNC) + cost_rad(2, e_2)) * 1.25
                print(f'Base cost (1.25 times the above cost):\t{RNC}')
                for e_3 in range(1, e_2 + 1, 1):
                    for e_5 in range(1, e_3 + 1, 1):
                        for e_7 in range(1, e_5 + 1, 1):
                            # ---
                            keyspace = KEYSPACE - float(security([e_2, e_3, e_5, e_7], 4))
                            e_tmp = list(e_init)[:(n - 3)]
                            assert(len(e_tmp) == (n - 3))
                            # Increasing to reach the right keyspace size
                            i = 0
                            while security(e_tmp, n - 3) < keyspace:
                                e_tmp[i] += 1
                                i = (i + 1) % (n - 3)
                            # Decreasing to reach the right keyspace size
                            i = 0
                            while security(e_tmp, n - 3) >= keyspace:
                                if e_tmp[i] > 0:
                                    e_tmp[i] -= 1
                                i = (i + 1) % (n - 3)

                            e_tmp = numpy.array(e_tmp + [0,0,0])
                            # ---
                            print("\nInitial integer vector of bounts (e_0, ..., e_%d)" % n)
                            printl("e", e_tmp, n // k)
                            RNC_prev, _, _, _, _ = strategy_block_cost(L[::-1], e_tmp)

                            print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RNC_prev[0] / (10.0**6), RNC_prev[1] / (10.0**6), RNC_prev[2] / (10.0**6), measure(RNC_prev) / (10.0**6)) )
                            print("\tSecurity ~ %f\n" % security(e_tmp, n))
                            print("_______________________________________________________________________________________________________________________________")
                            print("We proceed by searching a better integer vector of bounds\n")
                            for i in range(1, int(ceil( (1.0*setting.exponent) / (1.0*r) ))):
                                e_tmp, RNC_tmp = optimal_bounds(L[::-1], e_tmp, r, keyspace)
                                if measure(RNC_tmp) == measure(RNC_prev):
                                    break
                                else:
                                    RNC_prev = RNC_tmp

                            # Updating if the results were improved
                            RNC_tmp = measure(RNC_tmp) + cost_rad(2, e_2) + cost_rad(3, e_3) + cost_rad(5, e_5) + cost_rad(7, e_7)
                            if RNC_tmp >= RNC:
                                continue

                            assert(RNC_tmp < RNC)
                            RNC = RNC_tmp
                            e = e_tmp

                            print(f'Cost including the radical isogeny part:\t{RNC}')
                            print(f'Number of degree-2 isogeny constructions on the surface:\t{e_2}')
                            print(f'Number of degree-3 radical isogeny constructions:\t\t{e_3}')
                            print(f'Number of degree-5 radical isogeny constructions:\t\t{e_5}')
                            print(f'Number of degree-7 radical isogeny constructions:\t\t{e_7}')
                            # --------------------------------------------------------------------------------------------------
                            f = open("./tmp/" + setting.algorithm + "_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + raw + ".py", "w")
                            f.write( 'm = [' + ', '.join([ str(ei) for ei in e[::-1] ]) + ']')
                            if setting.algorithm == 'csurf':
                                f.write( '\ne_2 = %d' % e_2)
                                f.write( '\ne_3 = %d' % e_3)
                                f.write( '\ne_5 = %d' % e_5)
                                f.write( '\ne_7 = %d' % e_7)
                            f.close()
                            print("_______________________________________________________________________________________________________________________________\n")
            else:
                print(f'Cost assuming best CSIDH-configuration:\t{measure(RNC)}')
                RNC = measure(RNC) * 1.25
                print(f'Base cost (1.25 times the above cost):\t{RNC}')
                exp = 256
                assert(e_3 == 0)
                assert(e_5 == 0)
                assert(e_7 == 0)
                for e_2 in range(1, exp + 1, 1):
                    # ---
                    keyspace = KEYSPACE - float(security([e_2], 1))
                    e_tmp = list(e_init)
                    assert(len(e_tmp) == n)
                    i = 0
                    while security(e_tmp, n) >= keyspace:
                        if e_tmp[i] > 0:
                            e_tmp[i] -= 1
                        i = (i + 1) % n

                    e_tmp = numpy.array(e_tmp)
                    # ---
                    print("\nInitial integer vector of bounts (e_0, ..., e_%d)" % n)
                    printl("e", e_tmp, n // k)
                    RNC_prev, _, _, _, _ = strategy_block_cost(L[::-1], e_tmp)

                    print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RNC_prev[0] / (10.0**6), RNC_prev[1] / (10.0**6), RNC_prev[2] / (10.0**6), measure(RNC_prev) / (10.0**6)) )
                    print("\tSecurity ~ %f\n" % security(e_tmp, n))
                    print("_______________________________________________________________________________________________________________________________")
                    print("We proceed by searching a better integer vector of bounds\n")
                    for i in range(1, int(ceil( (1.0*setting.exponent) / (1.0*r) ))):
                        e_tmp, RNC_tmp = optimal_bounds(L[::-1], e_tmp, r, keyspace)
                        if measure(RNC_tmp) == measure(RNC_prev):
                            break
                        else:
                            RNC_prev = RNC_tmp

                    # Updating if the results were improved
                    RNC_tmp = measure(RNC_tmp) + cost_rad(2, e_2)
                    if RNC_tmp >= RNC:
                        continue

                    assert(RNC_tmp < RNC)
                    RNC = RNC_tmp
                    e = e_tmp

                    print(f'Cost including the degree-2 isogeny part:\t{RNC}')
                    print(f'Number of degree-2 isogeny constructions on the surface:\t{e_2}')
                    # --------------------------------------------------------------------------------------------------
                    f = open("./tmp/" + setting.algorithm + "_" + setting.prime  + "_" + setting.style + "_m" + str(setting.exponent) + raw + ".py", "w")
                    f.write( 'm = [' + ', '.join([ str(ei) for ei in e[::-1] ]) + ']')
                    if setting.algorithm == 'csurf':
                        f.write( '\ne_2 = %d' % e_2)
                        f.write( '\ne_3 = %d' % e_3)
                        f.write( '\ne_5 = %d' % e_5)
                        f.write( '\ne_7 = %d' % e_7)
                    f.close()
                    print("_______________________________________________________________________________________________________________________________\n")


if __name__ == "__main__":
    main()