# Don't forget the constant-time in CSURF!

Comparison between large constant-time implementations of CSIDH and CSURF using radical isogenies, in terms of field operations.


## Compilation

The syntax compilation can be viewed by running one of the following three commands

```bash
# Corresponding with the key-exchange protocol
python3 main.py -h or python3 main.py --help
# Corresponding with benchmarking (only for CSIDH, which has a variable running-time cost independent from the key)
python3 bench.py -h or python3 bench.py --help
```
and any of the above commands should show:

```bash
usage: main.py [-h] -a ALGORITHM -p PRIME -m MULTIEVAL [-s STYLE] [-b BENCHMARK] [-e EXPONENT] [-u] [-v] [-r]

Parses command.

optional arguments:
  -h, --help            show this help message and exit
  -a ALGORITHM, --algorithm ALGORITHM
                        csidh or csurf algorithm
  -p PRIME, --prime PRIME
                        prime number configuration should be stored in pSUFFIX (sop folder is taken as default).
  -m MULTIEVAL, --multieval MULTIEVAL
                        unscaled or scaled multi-evaluation procedure
  -s STYLE, --style STYLE
                        style to be used: wd1 (with dummy operations and a single torsion point), wd2 (with dummy operations and a two torsion point), or df (dummy-free approach (only for CSIDH)).
  -b BENCHMARK, --benchmark BENCHMARK
                        number of experiments to be used in the benchmark.
  -e EXPONENT, --exponent EXPONENT
                        For determining the number k of small odd primes to be used. The keyspace size is either (2e + 1)^k [wd2 style] or (e + 1)^n [wd1 and df styles].
  -u, --units           Used to precompute the running time of each small odd prime degree isogeny construction and evaluation with velusqrt formulae (required only at first run).
  -v, --verbose         Verbose mode.
  -r                    `raw' mode (for CSURF). With this argument, low degree isogenies are computed with the CSIDH group action evaluation. Without this, radical isogenies of these degrees are used.
```

Notice, `csidh` and `csurf`  option requires the use of the option `-s` (or `--style`). The verbose option `-v` (or `--verbose`) requires to have stored the tunned parameters `#I`, `#J` and `#K` in the new velu's formulae.

## Adding new primes

The field characteristic `p` should be stored in directory `sop/`.
Not all CSIDH primes are usable for CSURF: c should be 3 or larger (see below):

```bash
# CSIDH format (here p = 2^c * l_1 * .... l_n - 1)
c l_1 l_2 ... l_n
```
## Examples

We summarize some examples of runs as follows

```bash
# CSIDH + velusqrt + unscaled multievaluation approach + tunned velusqrt parameters
python3 main.py -p p512  -m unscaled -s wd2 -v -a csidh
python3 main.py -p p512  -m unscaled -s wd1 -v -a csidh
python3 main.py -p p512  -m unscaled -s df -v -a csidh
# CSIDH + velusqrt + scaled multievaluation approach + tunned velusqrt parameters
python3 main.py -p p512  -m scaled -s wd2 -v -a csidh
python3 main.py -p p512  -m scaled -s wd1 -v -a csidh
python3 main.py -p p512  -m scaled -s df -v -a csidh

# CSURF + velusqrt + unscaled multievaluation approach + tunned velusqrt parameters
python3 main.py -p p512  -m unscaled -s wd2 -v -a csurf -r
python3 main.py -p p512  -m unscaled -s wd1 -v -a csurf -r
python3 main.py -p p512  -m unscaled -s df -v -a csurf -r
# CSURF + velusqrt + scaled multievaluation approach + tunned velusqrt parameters
python3 main.py -p p512  -m scaled -s wd2 -v -a csurf -r
python3 main.py -p p512  -m scaled -s wd1 -v -a csurf -r
python3 main.py -p p512  -m scaled -s df -v -a csurf -r

# CRADS + velusqrt + unscaled multievaluation approach + tunned velusqrt parameters
python3 main.py -p p512  -m unscaled -s wd2 -v -a csurf
python3 main.py -p p512  -m unscaled -s wd1 -v -a csurf
python3 main.py -p p512  -m unscaled -s df -v -a csurf
# CRADS + velusqrt + scaled multievaluation approach + tunned velusqrt parameters
python3 main.py -p p512  -m scaled -s wd2 -v -a csurf
python3 main.py -p p512  -m scaled -s wd1 -v -a csurf
python3 main.py -p p512  -m scaled -s df -v -a csurf

# Benchmarking with 7 experiments
python3 bench.py -p p512  -m unscaled -s wd2 -v -a csidh -b 7
python3 bench.py -p p512  -m unscaled -s wd2 -v -a csurf -b 7
python3 bench.py -p p512  -m unscaled -s wd2 -v -a csurf -r -b 7
```

## Configuration

For the first instance execution, one must use the option `-u`, which precomputes the cost of the velusqrt formulae and then it computes the optimal strategies. Recall, the option `-v` is used to ensure the use of the tunned velusqrt parameters. However, If a new prime is going to be used, then one should run the following:

```bash
# Computing the tunned velusqrt parameters for the scaled multievaluation approach
bash tunned-parameters.sh p512 scaled
# Computing the tunned velusqrt parameters for the unscaled multievaluation approach
bash tunned-parameters.sh p512 unscaled

# Generating exponents bounds (e = 1), if e > 1 then one must add -e EXPONENT for some integer EXPONENT
# unscaled multievaluation
python3 bounds.py -p p512  -m unscaled -s wd2 -v -a csidh
python3 bounds.py -p p512  -m unscaled -s wd1 -v -a csidh
python3 bounds.py -p p512  -m unscaled -s df -v -a csidh

# scaled multievaluation
python3 bounds.py -p p512  -m scaled -s wd2 -v -a csidh
python3 bounds.py -p p512  -m scaled -s wd1 -v -a csidh
python3 bounds.py -p p512  -m scaled -s df -v -a csidh

# unscaled multievaluation (first run for generating all the auxiliar required files)
python3 main.py -p p512  -m unscaled -s wd2 -v -a csidh -u
python3 main.py -p p512  -m unscaled -s wd1 -v -a csidh -u
python3 main.py -p p512  -m unscaled -s df -v -a csidh -u
```

The optimal strategies are not changing for the unscaled and scaled approachs, for that reason we use the same strategies for boths. In other words, the difference in cost from the scaled and unscaled are negligeble for CSIDH and CSURF.

## Remarks

We focused on primes such that `p = 7 mod 8`. The author list is given in alphabetical ordering.


Per prime we use specific addition chains to compute the inverse and the N-th roots. These have been calculated using https://github.com/mmcloughlin/addchain
If you use different primes than the ones we are using, you need to convert those addition chains to be usable in the Python code.
The script for that is decoding_Addchains.py.

```bash
K=74
BITS=3072
python3 decoding_AddChains.py --bits $BITS --k $K --exponent inv
python3 decoding_AddChains.py --bits $BITS --k $K --exponent sq
python3 decoding_AddChains.py --bits $BITS --k $K --exponent tri
python3 decoding_AddChains.py --bits $BITS --k $K --exponent quart
python3 decoding_AddChains.py --bits $BITS --k $K --exponent quint
python3 decoding_AddChains.py --bits $BITS --k $K --exponent sept
python3 decoding_AddChains.py --bits $BITS --k $K --exponent novem
```

## Authors

1. **Jesús-Javier Chi-Domínguez** <jesus.chidominguez@tuni.fi>, <chidoys@gmail.com>; and
2. **Krijn Reijnders** <krijn.reijnders@ru.nl>, <reijnderskrijn@gmail.com>

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details.

## Funding

This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 804476). 
