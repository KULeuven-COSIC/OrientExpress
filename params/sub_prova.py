from sage.all import *

import sys 

if len(sys.argv) != 2:
    print("Usage: sage --python sub_prova.py a")
    exit()

a = int(sys.argv[1])

def is_nice(n):
    small_factors = []
    for c in range(2, 1000):
        while gcd(n, c) != 1:
            n //= gcd(n, c)
            small_factors.append(c)
    return is_prime(n), len(small_factors)

all_sols = []

for m in range(1, 100):
    if m%2 == 0:
        continue
    for n in range(m, 100):
        if n%2 == 0:
            continue
        QF = BinaryQF([n, 0, m])
        ls = QF.solve_integer(2**a, _flag=1)
        ls = [(abs(a), abs(b), m, n) for a,b in ls]
        ls = list(set(ls))
        all_sols += ls

for f1,f2,m,n in all_sols:
    is_nice_f1, n_small_fac_f1 = is_nice(f1)
    is_nice_f2, n_small_fac_f2 = is_nice(f2)
    n_fac_m = sum([x[1] for x in list(factor(m))])
    n_fac_n = sum([x[1] for x in list(factor(n))])
    if is_nice_f1 and is_nice_f2:
        print(f'{f1 = } {f2 = } {m = } {n = }')
        print(f'{factor(f1) = } {factor(f2) = }')
        print(f'used primes = {n_fac_m + n_fac_n + n_small_fac_f1 + n_small_fac_f2}')
        print()


# a = 71
# f1 = 10953205453 f2 = 2506940477 m = 13 n = 19 
# factor(f1) = 10953205453 factor(f2) = 139 * 18035543
# used primes = 3

# a = 101
# f1 = 161603664502629 f2 = 6642858156115 m = 47 n = 97
# factor(f1) = 3 * 23 * 2342082094241 factor(f2) = 5 * 1328571631223
# used primes = 5