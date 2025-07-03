from sage.all import *
proof.all(False)

import sys 

if len(sys.argv) != 5:
    print("Usage: sage --python sub_prova.py a f1 f2 n")
    exit()
# Find a, f1, f2 with sub_aa.py
a = int(sys.argv[1])
f1 = int(sys.argv[2])
f2 = int(sys.argv[3])
n = int(sys.argv[4]) # n = estimate number of "used" primes = #({divisors(m), divisors(n)}\cap{M_ls})

f = f1*f2 
assert f < 2**a
f_bits = ZZ(f).nbits()

if len(sys.argv) != 2:
    print("Usage: sage --python sub_prova.py n_primes")
    exit()
n_primes = ZZ(sys.argv[1])

r = n
M_bits_ls = []
for B in range(1,7):
    while ZZ((2*B+1)**(r-n)).nbits() < 256:
        r += 1
    print(f'{ZZ((2*B+1)**(r-n)).nbits() = }')
    M = primes_first_n(r)[1:]
    print(f'{M = }')
    last_prime = M[-1]
    M_bits = ceil(log(prod(M), 2))
    M_bits_ls.append(M_bits)
    cl_size = f_bits + (a + M_bits)//2
    print(f'{B = } {r = }')
    print(f'{last_prime = } {M_bits = }')
    print(f'{cl_size = } ')
    print()
    r = n
