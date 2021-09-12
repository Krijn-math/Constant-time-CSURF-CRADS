# Fully projective radical isogenies in constant-time

Comparison between large constant-time implementations of CSIDH and CSURF using (or not) radical isogenies, in terms of field operations.


## Dependencies and compiling the code

This Python-code implementation use the following python-package:

1. `progressbar`: `pip3 install progressbar`
1. `numpy`: `pip3 install numpy`
1. `functools` (it should be already installed from your python installation)
1. `argparse` (it should be already installed from your python installation)
1. `click`: `pip3 install click`

If you get into trouble when running the code, let us know.
The current version is working on Windows, Linux (Ubuntu-based systems), and macOS but handling all the systems is tricky.

The syntax compilation can be viewed by running one of the following three commands

```bash
# python3 main.py -h or python3 main.py --help
# python3 bench.py -h or python3 bench.py --help

$ python3 main.py --help
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
                        style to be used: wd1 (with dummy operations and a single torsion point), wd2 (with dummy operations and a two torsion point), or df (dummy-free approach).
  -b BENCHMARK, --benchmark BENCHMARK
                        number of experiments to be used in the benchmark.
  -e EXPONENT, --exponent EXPONENT
                        For determining the number k of small odd primes to be used. The keyspace size is either (2e + 1)^k [wd2 style] or (e + 1)^n [wd1 and df styles].
  -u, --units           Used to precompute the running time of each small odd prime degree isogeny construction and evaluation with velusqrt formulae (required only at first run).
  -v, --verbose         Verbose mode.
  -r, --radical         Radical isogenies mode (not only degree-2 isogenies on the surface.
```

Notice, the `--style df` (or `-s df`) option is not currently supported.
The verbose option `-v` (or `--verbose`) enables the use of square-root V&eacute;lu's tuned parameter.
Thus, it requires to have stored the tuned parameters `#I`, `#J` and `#K` for the square-root V&eacute;lu's formulas.

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
# Single random key-exchange

## CSIDH
### unscaled multievaluation approach
python3 main.py -p p512  -m unscaled -s wd2 -v -a csidh -e 5
python3 main.py -p p512  -m unscaled -s wd1 -v -a csidh -e 10
python3 main.py -p p1024  -m unscaled -s wd2 -v -a csidh -e 5
python3 main.py -p p1024  -m unscaled -s wd1 -v -a csidh -e 10
python3 main.py -p p1792  -m unscaled -s wd2 -v -a csidh -e 5
python3 main.py -p p1792  -m unscaled -s wd1 -v -a csidh -e 10
python3 main.py -p p2048  -m unscaled -s wd2 -v -a csidh -e 5
python3 main.py -p p2048  -m unscaled -s wd1 -v -a csidh -e 10
python3 main.py -p p3072  -m unscaled -s wd2 -v -a csidh -e 5
python3 main.py -p p3072  -m unscaled -s wd1 -v -a csidh -e 10
python3 main.py -p p4096  -m unscaled -s wd2 -v -a csidh -e 5
python3 main.py -p p4096  -m unscaled -s wd1 -v -a csidh -e 10
### scaled multievaluation approach
python3 main.py -p p512  -m scaled -s wd2 -v -a csidh -e 5
python3 main.py -p p512  -m scaled -s wd1 -v -a csidh -e 10
python3 main.py -p p1024  -m scaled -s wd2 -v -a csidh -e 5
python3 main.py -p p1024  -m scaled -s wd1 -v -a csidh -e 10
python3 main.py -p p1792  -m scaled -s wd2 -v -a csidh -e 5
python3 main.py -p p1792  -m scaled -s wd1 -v -a csidh -e 10
python3 main.py -p p2048  -m scaled -s wd2 -v -a csidh -e 5
python3 main.py -p p2048  -m scaled -s wd1 -v -a csidh -e 10
python3 main.py -p p3072  -m scaled -s wd2 -v -a csidh -e 5
python3 main.py -p p3072  -m scaled -s wd1 -v -a csidh -e 10
python3 main.py -p p4096  -m scaled -s wd2 -v -a csidh -e 5
python3 main.py -p p4096  -m scaled -s wd1 -v -a csidh -e 10

## CSURF
### unscaled multievaluation approach
python3 main.py -p p512  -m unscaled -s wd2 -v -a csurf -e 5
python3 main.py -p p512  -m unscaled -s wd1 -v -a csurf -e 10
python3 main.py -p p1024  -m unscaled -s wd2 -v -a csurf -e 5
python3 main.py -p p1024  -m unscaled -s wd1 -v -a csurf -e 10
python3 main.py -p p1792  -m unscaled -s wd2 -v -a csurf -e 5
python3 main.py -p p1792  -m unscaled -s wd1 -v -a csurf -e 10
python3 main.py -p p2048  -m unscaled -s wd2 -v -a csurf -e 5
python3 main.py -p p2048  -m unscaled -s wd1 -v -a csurf -e 10
python3 main.py -p p3072  -m unscaled -s wd2 -v -a csurf -e 5
python3 main.py -p p3072  -m unscaled -s wd1 -v -a csurf -e 10
python3 main.py -p p4096  -m unscaled -s wd2 -v -a csurf -e 5
python3 main.py -p p4096  -m unscaled -s wd1 -v -a csurf -e 10
### scaled multievaluation approach
python3 main.py -p p512  -m scaled -s wd2 -v -a csurf -e 5
python3 main.py -p p512  -m scaled -s wd1 -v -a csurf -e 10
python3 main.py -p p1024  -m scaled -s wd2 -v -a csurf -e 5
python3 main.py -p p1024  -m scaled -s wd1 -v -a csurf -e 10
python3 main.py -p p1792  -m scaled -s wd2 -v -a csurf -e 5
python3 main.py -p p1792  -m scaled -s wd1 -v -a csurf -e 10
python3 main.py -p p2048  -m scaled -s wd2 -v -a csurf -e 5
python3 main.py -p p2048  -m scaled -s wd1 -v -a csurf -e 10
python3 main.py -p p3072  -m scaled -s wd2 -v -a csurf -e 5
python3 main.py -p p3072  -m scaled -s wd1 -v -a csurf -e 10
python3 main.py -p p4096  -m scaled -s wd2 -v -a csurf -e 5
python3 main.py -p p4096  -m scaled -s wd1 -v -a csurf -e 10

# CRADS (CSURF  + radical isogenies)
### unscaled multievaluation approach
python3 main.py -p p512  -m unscaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p512  -m unscaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p1024  -m unscaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p1024  -m unscaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p1792  -m unscaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p1792  -m unscaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p2048  -m unscaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p2048  -m unscaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p3072  -m unscaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p3072  -m unscaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p4096  -m unscaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p4096  -m unscaled -s wd1 -v -a csurf -e 10 -r
### scaled multievaluation approach
python3 main.py -p p512  -m scaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p512  -m scaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p1024  -m scaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p1024  -m scaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p1792  -m scaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p1792  -m scaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p2048  -m scaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p2048  -m scaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p3072  -m scaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p3072  -m scaled -s wd1 -v -a csurf -e 10 -r
python3 main.py -p p4096  -m scaled -s wd2 -v -a csurf -e 5 -r
python3 main.py -p p4096  -m scaled -s wd1 -v -a csurf -e 10 -r

# Benchmarking with 7 experiments
python3 bench.py -p p512  -m unscaled -s wd2 -v -a csidh -b 7 -e 5
python3 bench.py -p p512  -m unscaled -s wd2 -v -a csurf -b 7 -e 5
python3 bench.py -p p512  -m unscaled -s wd2 -v -a csurf -r -b 7 -e 5
```

## Configuration

For the first instance execution, one must use the option `-u`, which precomputes the cost of the square-root V&eacute;lu formulas, and then it computes the optimal strategies.
Recall, the option `-v` is used to ensure the use of the tuned square-root V&eacute;lu parameters.
However, If a new prime is going to be used, then one should run the following:

```bash
# Computing the tunned velusqrt parameters for the scaled multievaluation approach
bash tunned-parameters.sh p512 scaled
# Computing the tunned velusqrt parameters for the unscaled multievaluation approach
bash tunned-parameters.sh p512 unscaled

# Generating exponents bounds (next runs can takes hours, even days depending on the search space)
# unscaled multievaluation
python3 bounds.py -p p512  -m unscaled -s wd2 -v -a csidh -e 5 
python3 bounds.py -p p512  -m unscaled -s wd1 -v -a csidh -e 10

# scaled multievaluation
python3 bounds.py -p p512  -m scaled -s wd2 -v -a csidh
python3 bounds.py -p p512  -m scaled -s wd1 -v -a csidh

# Example: CRADS + unscaled multievaluation (first run for generating all the auxiliar required files)
ALGORITHM=csurf
CRADS=-r

python3 main.py -p p512  -m unscaled -s wd2 -v -a $ALGORITHM -e 5 -u $CRADS
python3 main.py -p p512  -m unscaled -s wd1 -v -a $ALGORITHM -e 10 -u $CRADS
python3 main.py -p p1024  -m unscaled -s wd2 -v -a $ALGORITHM -e 5 -u $CRADS
python3 main.py -p p1024  -m unscaled -s wd1 -v -a $ALGORITHM -e 10 -u $CRADS
python3 main.py -p p1792  -m unscaled -s wd2 -v -a $ALGORITHM -e 5 -u $CRADS
python3 main.py -p p1792  -m unscaled -s wd1 -v -a $ALGORITHM -e 10 -u $CRADS
python3 main.py -p p2048  -m unscaled -s wd2 -v -a $ALGORITHM -e 5 -u $CRADS
python3 main.py -p p2048  -m unscaled -s wd1 -v -a $ALGORITHM -e 10 -u $CRADS
python3 main.py -p p3072  -m unscaled -s wd2 -v -a $ALGORITHM -e 5 -u $CRADS
python3 main.py -p p3072  -m unscaled -s wd1 -v -a $ALGORITHM -e 10 -u $CRADS
python3 main.py -p p4096  -m unscaled -s wd2 -v -a $ALGORITHM -e 5 -u $CRADS
python3 main.py -p p4096  -m unscaled -s wd1 -v -a $ALGORITHM -e 10 -u $CRADS
```

The optimal strategies are not changing for the unscaled and scaled approaches, for that reason we use the same strategies for boths.
In other words, the difference in cost from the scaled and unscaled are negligible for CSIDH and CSURF.

## Remarks

We focused on primes such that `p = 7 mod 8`. The author list is given in alphabetical ordering.

Per prime we use specific addition chains to compute the inverse, and the N-th roots.
These have been calculated using https://github.com/mmcloughlin/addchain
If you use different primes than the ones we are using, you need to convert those addition chains to be usable in the Python-code.
The script for that is `decoding_Addchains.py` (it use `click` package).

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

<!--
## Authors

1. **Jesús-Javier Chi-Domínguez** <jesus.dominguez@tii.ae>, <chidoys@gmail.com>; and
2. **Krijn Reijnders** <krijn.reijnders@ru.nl>, <reijnderskrijn@gmail.com>
-->

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details.

## Funding

This project has initially received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 804476). 
