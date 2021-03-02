from framework import *

try:
    exec("from %s import *" % setting.algorithm)

except:

    print("  ,-~~-.___.          ")
    print(" / |  '     \\          SYNTAX ERROR ..., run python3 %s -h for help" % sys.argv[0]) 
    print("(  )         0        ")  
    print(" \_/-, ,----'         ")         
    print("    ====           // ")
    print("   /  \-'~;    /~~~(O)")
    print("  /  __/~|   /       |")   
    print("=(  _____| (_________|")
    exit(7)

print("// The exponents bounds and private keys go from the smallest to the largest small odd primes ell_i\'s\n")
if (setting.algorithm == 'csurf'):
    if not setting.raw:
        privkey_print('Exponents bounds', [e_2] + [e_3] + [e_5] + [e_7] + m[3:])
    else:
        privkey_print('Exponents bounds', [e_2] + m)
else:
    privkey_print('Exponents bounds', m)

print("")
tmp = list(set(m))
SAMPLE = [ [0.0, 0.0, 0.0] ] * setting.benchmark
SAMPLE_RAD = [ [0.0, 0.0, 0.0] ] * setting.benchmark
SAMPLE_VELU = [ [0.0, 0.0, 0.0] ] * setting.benchmark
SAMPLE_VALIDATE = [ [0.0, 0.0, 0.0] ] * setting.benchmark
bar = Bar('// Running experiments' , max=setting.benchmark)

A = 0   # E : y² = x³ + x
for main_i in range(setting.benchmark):

    if (setting.algorithm == 'csurf'):
        if not setting.raw:
            a = [ random.randint(0, e_i) for e_i in [e_2, e_3, e_5, e_7] ] + random_key(m)[3:] # Random private key
        else:
            a = [random.randint(0, e_2)] + random_key(m) # Random private key
    else:
        a = random_key(m)                       # Random private key
    A, SAMPLE_VELU[main_i], SAMPLE_RAD[main_i], SAMPLE[main_i] = derive(a, A)        # Random public curve E : y² = x³ + Ax² + x
    set_zero_ops()
    V = validate(A)
    assert(V)
    SAMPLE_VALIDATE[main_i] = get_ops()
        
    #print("Random(EllipticCurve(x^3 + 0x%X * x^2 + x)) * (p+1);" % coeff(B))
        
    bar.next()
bar.finish()

AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE ]), statistics.mean([ ops[1] for ops in SAMPLE ]), statistics.mean([ ops[2] for ops in SAMPLE ]) ]
AVERAGE_RAD = [statistics.mean([ ops[0] for ops in SAMPLE_RAD ]), statistics.mean([ ops[1] for ops in SAMPLE_RAD ]), statistics.mean([ ops[2] for ops in SAMPLE_RAD ]) ]
AVERAGE_VELU = [statistics.mean([ ops[0] for ops in SAMPLE_VELU ]), statistics.mean([ ops[1] for ops in SAMPLE_VELU ]), statistics.mean([ ops[2] for ops in SAMPLE_VELU ]) ]
print("\n// Average number of field operations (RAD):\t\t%2.3fM + %2.3fS + %2.3fa := %2.3fM" % (AVERAGE_RAD[0] / (10.0**6), AVERAGE_RAD[1] / (10.0**6), AVERAGE_RAD[2] / (10.0**6), measure(AVERAGE_RAD) / (10.0**6)) )
print("// Average number of field operations (GAE):\t\t%2.3fM + %2.3fS + %2.3fa := %2.3fM" % (AVERAGE_VELU[0] / (10.0**6), AVERAGE_VELU[1] / (10.0**6), AVERAGE_VELU[2] / (10.0**6), measure(AVERAGE_VELU) / (10.0**6)) )
print("// Average number of field operations (total):\t\t%2.3fM + %2.3fS + %2.3fa := %2.3fM" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )

AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE_VALIDATE ]), statistics.mean([ ops[1] for ops in SAMPLE_VALIDATE ]), statistics.mean([ ops[2] for ops in SAMPLE_VALIDATE ]) ]
print("// Average number of field operations (validate):\t%2.3fM + %2.3fS + %2.3fa := %2.3fM\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )


