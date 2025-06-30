import sys
sys.path.insert(0, '..')

from sage.all import *

from isolib.utils import BiDLP
from isolib.utils import weil_pairing_pari
from isolib.utils import torsion_basis
from isolib.utils import FactIso
from isolib.HDiso import HDChain

def smooth_factor(n, B=1000):
    """
    Remove factors up to a bound and check the rest for primality.
    """
    for ell in primes(2, B):
        v_l = valuation(n, ell)
        if v_l > 0:
            n //= ell**v_l

    if is_prime(n):
        return n
    return False

def RepresentInteger(M, O0):
    """
    Returns random quaternions of reduced norm M in O0
    """
    B = O0.quaternion_algebra()
    i, j, k = B.gens()

    p = ZZ(-j**2) # p
    q = ZZ(-i**2) # 1

    QF = BinaryQF([1, 0, q]) # 1*x^2 + 0*xy + q*y^2

    # Buond from SQISign
    m = isqrt(M/(p*(1+q)))

    n_att = m**3 # Some margin
    # print(f'{m = }')
    for _ in range(n_att):
        z = randint(1, m)
        t = randint(1, m)

        Mm = M - p*QF(z, t) # QF(z, t) == (z + i*t).reduced_norm()
        ell = smooth_factor(Mm)
        if not ell:
            continue
        pari.addprimes(ell)
        sol = QF.solve_integer(Mm) # TODO: Mm should be easy to factor
        pari.removeprimes(ell)
        if sol:
            return sol[0] + i*sol[1] + j*(z + i*t)

    return None

def ideal_to_kernel_E0(theta, l, e, D):
    """
    Given a quaternion element theta, find a generator of its kernel on the l^e
    torsion. D contains the curve E0 and maps corresponding to i and j.
    """
    print(f'Finding kernel for {l = } {e = }')
    assert theta.denominator() == 1, "fractional theta"
    E0, qt_i, qt_j = D

    # Evaluate theta on both points
    P, Q, ePQ = torsion_basis(E0, l, e)
    print('\t- Torsion basis done')
    thetaP = eval_end(theta, P, D)
    thetaQ = eval_end(theta, Q, D)
    print('\t- Point evaluation done')

    # Find kernel matrix and points
    a, b = BiDLP(thetaP, P, Q, l**e, ePQ)
    c, d = BiDLP(thetaQ, P, Q, l**e, ePQ)

    if e == 1:
        M = Matrix(Integers(l**e), 2, [a, b, c, d])
        s, t = M.left_kernel().gen()
    else:
        # Pallepiene version of echelon form - fails if all the matrix is 0 mod l
        if int(c) % l != 0 or int(a) % l != 0:
            s, t = c, -a
        # elif int(d) % l != 0 or int(b) % l != 0:
        else:
            s, t = d, -b
    K = s*P + t*Q
    print('\t- Matrix done')

    # Sanity checks
    assert eval_end(theta, K, D) == 0, "point not in kernel"
    assert (l**e)*K == 0, "point of wrong order"
    # WARNING: K may not be full order (especially for 2)
    return K

def eval_end(theta, P, D, valP=None):
    """
    Evaluate the endomorpshism theta in End(E0) on P. D contains E0 and the
    maps corresponding to i and j. Optionally, valP contains the precomputation
    of values [P, iP, jP, kP].
    """
    if P == 0:
        return P
    E0, qt_i, qt_j = D

    if not valP:
        iP = qt_i(P)
        jP = qt_j(P)
        kP = qt_i(jP)
        valP = [P, iP, jP, kP]

    ell = P.order()
    # print(f'{factor(P.order()) = }')
    thetaP = sum([(ZZ(c) % ell)*d for c, d in zip(list(theta), valP)])
    return thetaP

# B = 2 
# a = 1184 
# f1 = 1231254475685351434422380365578760673087925822404846536968852569977683419024788296803556490898281430118768907794888389473420330169144416825302777123752230103113816123657832222249 
# f2 = 75033646116264288529343216199419546564688109171625450229668622752873962650479131951565410744026425835324324084605894029954440843233851790981334170547193707616293911428648033111 
# m = 173 
# n = 83 
# n_primes = 97

# B = 5
# a = 333 
# f1 = 11316416821907091267376074948469150738154508493417 
# f2 = 74961544251765932689098739581437248473675373596757 
# m = 5 
# n = 3
# n_primes = 66

# for n_primes use output "r" of sub_p.py


if len(sys.argv) != 8:
    print("Usage: sage --python sub_pgen.py B a f1 f2 m n n_primes")
    exit()

B = int(sys.argv[1])
a = int(sys.argv[2])
f1 = int(sys.argv[3])
f2 = int(sys.argv[4])
m = int(sys.argv[5]) 
n = int(sys.argv[6]) 
n_primes = int(sys.argv[7]) 
assert m*f1**2 + n*f2**2 == 2**a

f = f1*f2

M_ls = primes_first_n(n_primes+1)[1:]
print(f'{M_ls = }')
M = product(M_ls) 

# Find p
cof = 1
div = M*2**a
while True:
    cof += 2 # Odd cofactor so the T-isogeny does not interfere with 2-torsion
    if gcd(m*n, cof) == 1 and is_prime(cof*div - 1):
        break 
p = cof*div - 1
print(f'Found {p = }')
print(f'{cof = }')

# Order evaluation
eval_even = a
eval_odd = M // (m*n)

# Compute eigenvalues 
eigenvals = {}
for mi in M_ls:
    if f%mi == 0:
        continue
    R, x = PolynomialRing(GF(mi), names='x').objgen()
    rts = (x**2 + (f**2)*p).roots(multiplicities=False)
    # print(f'{mi = } {rts = }')
    assert len(rts) == 2 and rts[0] == -rts[1] # sanity check
    eigenvals[int(mi)] = [int(zz) for zz in rts]

# Find starting curve and tau by finding an endomorphism theta of E0 
# factoring through an f-isogeny 

g = 2**a - f
assert g > 0

E = EllipticCurve(Fp2, [1,0])

from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
quat_i = WeierstrassIsomorphism(E, (ii, 0, 0, 0)) # (x, y) -> (-x, iy)
quat_j = E.frobenius_isogeny()

Quat,(i,j,k) = QuaternionAlgebra(-1, -p).objgens()
O0 = Quat.quaternion_order([1, i, (i+j)/2, (1+k)/2])

D = (E, quat_i, quat_j)
P0, Q0, _ = torsion_basis(E, 2, a)

if f*g <= p:
    print(f'RepresentInteger with auxiliary primes')
    # RepresentIntegern with norm f*(2**a - f)*primes
    T = cof*M
    while True:
        theta = RepresentInteger(f*g*T, O0) # Need to make sure degree is big enough # Use improved RepresentInteger
        print(f'{theta = }')
        # if not theta:
        #     raise ValueError(f'No output from RepresentInteger')
        theta_bar = theta.conjugate()
        assert theta.reduced_norm() == f*g*T
        print('Finding kernel')
        ells = list(factor(T))
        K = E(0)
        fail = False
        for l, e in ells:
            Kl = ideal_to_kernel_E0(theta_bar, l, e, D)
            Kl.set_order(multiple = l**(e+1))
            ordK = Kl.order()

            if ordK != l**e:
                print(f'Failed for {l = }')
                fail = True # theta factors through mult by l (probably)
                break
            K += Kl
        if fail:
            continue

        assert eval_end(theta_bar, K, D) == 0
        break
    print({f'{K = }'})

    # Evaluate the T-isogeny corresponding to K
    print('Computing T-isogeny')
    phiT_hat = FactIso(K, d=T, fac=ells) # Warning: may take long
    assert phiT_hat(K) == 0
    E1 = phiT_hat.codomain()

    # Remove the T part from the endomorphism
    P = T*P0 # T here comes from applying the dual
    Q = T*Q0
    thetaP = eval_end(theta, P0, D)
    thetaQ = eval_end(theta, Q0, D)
    phidP = phiT_hat(thetaP)  # phid is the part of theta of deg = f*(2**a-f)
    phidQ = phiT_hat(thetaQ) 

    # The 2D isogeny Phi of kernel ((f*P, phidP), (f*Q, phidQ)) has our wanted psi of degree f as a component
    # Locating the factor of degree f in Phi
    print(f'Computing the chain')
    K1 = (f*P, phidP)
    K2 = (f*Q, phidQ)
    print(f'{K1[0].weil_pairing(K2[0],2**a)*K1[1].weil_pairing(K2[1],2**a) = }')
    Phi_chain = HDChain((K1,K2), a)
    P1, Q1 = Phi_chain.eval_two_points((P,E1(0)),(Q,E1(0)))
    if P.weil_pairing(Q, Q.order())**f == P1[1].weil_pairing(Q1[1], Q1[1].order()):
        # print(f"f-isogeny in component 2,1")
        x = 2
    elif P.weil_pairing(Q, Q.order())**f == P1[0].weil_pairing(Q1[0], Q1[0].order()):
        # print(f"f-isogeny in component 1,1")
        x = 1
    else:
        raise ValueError(f"weil pairing on image points of Phi_chain does not match")

else:
    print(f'RepresentInteger without auxiliary')
    theta = RepresentInteger(f*g, O0)
    print(f'{theta = }')
    P = P0
    Q = Q0
    thetaP = eval_end(theta, P0, D)
    thetaQ = eval_end(theta, Q0, D)
    phidP = thetaP
    phidQ = thetaQ

    # The 2D isogeny Phi of kernel ((f*P, phidP), (f*Q, phidQ)) has our wanted psi of degree f as a component
    # Locating the factor of degree f in Phi
    print(f'Computing the chain')
    K1 = (f*P, phidP)
    K2 = (f*Q, phidQ)
    print(f'{K1[0].weil_pairing(K2[0],2**a)*K1[1].weil_pairing(K2[1],2**a) = }')
    Phi_chain = HDChain((K1,K2), a)
    P1, Q1 = Phi_chain.eval_two_points((P,E(0)),(Q,E(0)))
    if P.weil_pairing(Q, Q.order())**f == P1[1].weil_pairing(Q1[1], Q1[1].order()):
        # print(f"f-isogeny in component 2,1")
        x = 2
    elif P.weil_pairing(Q, Q.order())**f == P1[0].weil_pairing(Q1[0], Q1[0].order()):
        # print(f"f-isogeny in component 1,1")
        x = 1
    else:
        raise ValueError(f"weil pairing on image points of Phi_chain does not match")

# Sanity
print('Sanity')
assert P.order() == 2**a
assert Q.order() == 2**a
assert P1[x-1].order() == 2**a
assert Q1[x-1].order() == 2**a

E = Phi_chain.codomain()[x-1] # codomain should be a product isogeny
# E_j_invariant = [int(i) for i in E.j_invariant().list()]
E_coeffs = [[int(c[0]), int(c[1])] for c in [E.a1(),E.a2(),E.a3(),E.a4(),E.a6()]]

# P0, Q0, _ = torsion_basis(E0.E, 2, a) # better than E0.torsion_basis(2**a)
# K1 = (f*P0, eval_end(theta, P0, D))
# K2 = (f*Q0, eval_end(theta, Q0, D))
# Phi_chain = HDChain((K1,K2), a)

# # The 2D isogeny Phi of kernel (f*P0, f*Q0), (theta(P0), theta(Q0)) has our wanted psi of degree f as a component
# # Locating the factor of degree f in Phi
# P1, Q1 = Phi_chain.eval_two_points((P0,0),(Q0,0))
# if P0.weil_pairing(Q0, Q0.order())**f == P1[1].weil_pairing(Q1[1], Q1[1].order()):
#     # print(f"f-isogeny in component 2,1")
#     x = 2  
# elif P0.weil_pairing(Q0, Q0.order())**f == P1[0].weil_pairing(Q1[0], Q1[0].order()):
#     # print(f"f-isogeny in component 1,1")
#     x = 1
# else:
#     raise ValueError(f"weil pairing on image points of Phi_chain does not match")

# E = Phi_chain.codomain()[x-1] # codomain should be a product isogeny
# # E_j_invariant = [int(i) for i in E.j_invariant().list()]
# E_coeffs = [[int(c[0]), int(c[1])] for c in [E.a1(),E.a2(),E.a3(),E.a4(),E.a6()]]


# Image of 2**a-torsion via psi
psi_P = P1[x-1] 
psi_Q = Q1[x-1]
print(f'P.weil_pairing(Q, 2**a)**f == psi_P.weil_pairing(psi_Q, 2**a) = {P.weil_pairing(Q, 2**a)**f == psi_P.weil_pairing(psi_Q, 2**a)}')
 
# points for test at the end "theta_j does not generate -new- isogeny"
psi_P_conj = E.frobenius_isogeny()(psi_P)
psi_Q_conj = E.frobenius_isogeny()(psi_Q)

# Dual of psi
print(f'Computing dual of psi')
U, V, _ = torsion_basis(E, 2, a)
u1, u2 = BiDLP(U, psi_P, psi_Q, 2**a)
v1, v2 = BiDLP(V, psi_P, psi_Q, 2**a)
psi_dual_U = u1*f*P + u2*f*Q
psi_dual_V = v1*f*P + v2*f*Q
print(f'U.weil_pairing(V, 2**a)**f == psi_dual_U.weil_pairing(psi_dual_V, 2**a) = {U.weil_pairing(V, 2**a)**f == psi_dual_U.weil_pairing(psi_dual_V, 2**a)}')
# print(f'log(psi_dual_U.weil_pairing(psi_dual_V, 2**a), U.weil_pairing(V, 2**a)) = {log(psi_dual_U.weil_pairing(psi_dual_V, 2**a), U.weil_pairing(V, 2**a))}')

# Computing the the frobenius conjugate of psi composed with psi dual on U,V
print(f'Computing tau')
R = quat_j(psi_dual_U) # apply frobenius to psi_dual_U
S = quat_j(psi_dual_V) 
print(f'quat_j(R) == psi_dual(U)  = {quat_j(R) == psi_dual_U}')
print(f'quat_j(S) == psi_dual(V) = {quat_j(S) == psi_dual_V}')
r1, r2 = BiDLP(R, P, Q, 2**a)
s1, s2 = BiDLP(S, P, Q, 2**a)
psi_R = r1*psi_P + r2*psi_Q
psi_S = s1*psi_P + s2*psi_Q
# print(f'{psi_R = } {psi_S = }')
E_frob = E.frobenius_isogeny()
img_U = E_frob(psi_R) # note frob(psi(R)) = psi_conj(frob(R)) = psi_conj(psi_dual_U)
img_V = E_frob(psi_S)

print(f'U.weil_pairing(V, 2**a)**(f**2) == img_U.weil_pairing(img_V, 2**a) = {U.weil_pairing(V, 2**a)**(f**2) == img_U.weil_pairing(img_V, 2**a)}')
print(f'log(U.weil_pairing(V, 2**a), img_U.weil_pairing(img_V, 2**a)**(f**2)) == -1 = {U.weil_pairing(V, 2**a)**(-f**2) == img_U.weil_pairing(img_V, 2**a)}')
if img_U == u1*f*psi_P_conj + u2*f*psi_Q_conj or img_V == v1*f*psi_P_conj + v2*f*psi_Q_conj:
    raise ValueError(f"theta_j does not generate -new- isogeny")

tau = [int(i) for i in img_U.x().list() + img_U.y().list() + img_V.x().list() + img_V.y().list()]
E_tors = [int(i) for i in U.x().list() + U.y().list() + V.x().list() + V.y().list()]

# print(f'{f1 = } {f2 = } {f = }')
# print(f'primes = {M_ls}')
# print(f'{p = }')
# print(f'{cof = }')
# print(f'{eval_even = }')
# print(f'{factor(eval_odd) = }')
# print(f'{m = }')
# print(f'{n = }')
# print(f'{E = } ')
# print(f'{E_tors = }')
# print(f'{tau = }')
# print(f'{eigenvals = }')

param_set = {
        'a':int(a),
        'B':int(B),
        'p':int(p), 
        'cof':int(cof),
        'eval_even':int(eval_even),
        'eval_odd':int(eval_odd),
        'f1':int(f1), 'f2':int(f2), 'f':int(f),
        'm':int(m), 'n':int(n),
        'E':E_coeffs, 
        'E_tors':E_tors,
        'tau':tau, 
        'primes':[int(i) for i in M_ls],
        'eigenvals':eigenvals
        }

from json import dump
dump(param_set, open(f'sub_pgen{a}.json', 'w'))

def check_dump(param_set):
    p = param_set['p']
    Fp2, (ii, ) = GF((p, 2), names='ii', modulus=[1,0,1]).objgens()
    E_coeffs = [param_set['E'][c][0] + ii*param_set['E'][c][1] for c in range(0,len(param_set['E']))]
    # E_coeffs = [param_set['E'][c][0] + ii*param_set['E'][c][1] for c in range(0,5)]
    E_dump = EllipticCurve(Fp2, E_coeffs)
    test1 = E_dump == E
    # print(f"{E_dump = }")
    # print(f"{E = }")
    print(f"{test1 = }")

    tau_U_V = [param_set['tau'][2*c] + ii*param_set['tau'][2*c+1] for c in range(0,len(param_set['tau'])//2)]
    E_conj = E_dump.frobenius_isogeny().codomain()
    points = [E_conj(tau_U_V[2*d], tau_U_V[2*d+1]) for d in range(0,len(tau_U_V)//2)]
    point_1 = points[0]
    point_2 = points[1]
    test2 = (point_1 == img_U and point_2 == img_V)
    # print(f"{point_1 = }")
    # print(f"{img_U = }")
    print(f"{test2 = }")

    E_dump_tors = [param_set['E_tors'][2*c] + ii*param_set['E_tors'][2*c+1] for c in range(0,len(param_set['E_tors'])//2)]
    E_U_V = [E_dump(E_dump_tors[2*d], E_dump_tors[2*d+1]) for d in range(0, len(E_dump_tors)//2)]
    U_dump = E_U_V[0]
    V_dump = E_U_V[1]
    test3 = U_dump.weil_pairing(V_dump, 2**(param_set['a']))**(f**2) == point_1.weil_pairing(point_2, 2**(param_set['a']))
    print(f"{test3 = }")

    return (test1 and test2) and test3

assert check_dump(param_set)
print(E_coeffs)

    


