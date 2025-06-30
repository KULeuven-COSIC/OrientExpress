import sys

from sage.all import *
proof.all(False)

# 1500 set
B = 5
n_primes = 64
a = 360
max_missing_primes = 15

# 4000 set
B = 2
a = 1102
n_primes = 95
max_missing_primes = 2

assert a % 2 == 0
M_ls = primes_first_n(n_primes+1)[1:]
M = product(M_ls)

# Find p
cof = 1
div = M * 2**(a+1)
while not is_prime(cof*div - 1):
    print(f'- {cof = }')
    cof += 1

p = cof*div - 1
print(f'Found {p = }')
print(f'-> {cof = }')


# Maximum number of primes that we do not include in x here
for miss in range(max_missing_primes):
    if miss == 0:
        used_primes = M_ls[:]
        other_primes = []
    else:
        used_primes = M_ls[:-miss]
        other_primes = M_ls[-miss:]
    print(f'Trying with {miss = } -> {other_primes}')
    Mx = product(used_primes)
    print(f'{Mx.nbits() = } {a = }')
    prime_pairs_found = 0

    # Find d
    l1 = ZZ(2**a)
    l2 = ZZ(2**a)
    x0 = 0

    d = None
    while True:
        if l1 < Mx:
            break
        x0 += 1
        if x0 % 10000 == 0:
            print(f'{x0 = }')
        # print(f'- {x0 = }')
        l1 -= Mx
        l2 += Mx

        if is_prime(l1) and is_prime(l2):
            # Check if the remeaning primes are split
            print('Prime pair found')
            prime_pairs_found += 1
            print(f'Found: {prime_pairs_found} - exp: {2**miss}')
            d1 = l1*l2
            is_ok = True
            for mi in other_primes:
                if kronecker(d1, mi) != 1:
                    is_ok = False
                    break
            if is_ok:
                print('Full pair found')
                print(f'{x0 = }')
                print(f'{Mx = }')
                print(f'{miss = }')
                d = d1
                break
            else:
                print('Fail')

    if d:
        break

if not d:
    print('AAAA')
    exit()

# Set convention on eigenvalues
eigenvals = {}
for mi in M_ls:
    R, x = PolynomialRing(GF(mi), names='x').objgen()
    rts = (x**2 + l1*l2*p).roots(multiplicities=False)
    assert len(rts) == 2 and rts[0] == -rts[1] # sanity check
    eigenvals[int(mi)] = [int(zz) for zz in rts]

# Find alpha and beta such that d = norm(alpha + beta*omega)
pari.addprimes([ZZ(l1), ZZ(l2)])
Q = BinaryQF([1, -1, 1])
alpha, beta = Q.solve_integer(l1*l2)

# Optionally: find elements of norm l1, l2
# alpha1, beta1 = Q.solve_integer(l1)
# alpha2, beta2 = Q.solve_integer(l2)
# print(f'{alpha1 = } {beta1 = }')
# print(f'{beta1 % 3 = } {beta1 % 5 = } {beta1 % 7 = }')
# print(f'{alpha2 = } {beta2 = } {beta2}')
# print(f'{beta2 % 3 = } {beta2 % 5 = } {beta2 % 7 = }')

print(f'{l1 = } {l2 = }')
print(f'primes = {M_ls}')
print(f'{p = }')
print(f'{alpha = } {beta = }')
print(f'{cof = }')
print(f'{x0 = }')
print(f'Last prime: {M_ls[-1]}')
#  print(f'{eigenvals = }')
h = (d*p).nbits()
print(f'Prime size: {p.nbits()} bits')
print(f'(Square) class group size: {h} bits')

param_set = {
        'a':int(a),
        'B':int(B),
        'p':int(p), 'cof':int(cof),
        'l1':int(l1), 'l2':int(l2),
        'alpha':int(alpha), 'beta':int(beta),
        'primes':[int(i) for i in M_ls],
        'eigenvals':eigenvals
        }

from json import dump
dump(param_set, open(f'cs_ftf{h}.json', 'w'))

