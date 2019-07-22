import argparse
import itertools

parser = argparse.ArgumentParser(description='Check hyperbolic isometries below an affine Coxeter element.')
parser.add_argument('case', type=str, help='one of F4, E6, E7, E8')
parser.add_argument('-v', '--verbose', action='store_true')

args = parser.parse_args()
print "Case {}".format(args.case)

family = args.case[0]
assert family in ['C', 'E', 'F']
m = int(args.case[1])

if family == 'E':
    assert m in [6, 7, 8]
    n = 8   # dimension of the ambient space
elif family == 'F':
    assert m == 4
    n = 4   # dimension of the ambient space
else:
    n = m   # dimension of the ambient space


# ambient space
V = VectorSpace(QQ, n)
Aff = AffineGroup(n, QQ)
# id = Isometry(Aff.one(), length=0)

e = V.basis()
assert all(e[i] == vector(matrix.identity(n)[i,:]) for i in xrange(n))


# root systems [Hum92, Section 2.10]
# (we only list one representative for each pair of opposite roots)
if family == 'E':
    roots = [e[i] - e[j] for i in xrange(n) for j in xrange(i+1, n)] + \
        [e[i] + e[j] for i in xrange(n) for j in xrange(i+1, n)]

    for comb in itertools.product([1,-1], repeat=n-1):
        if comb.count(-1) % 2 == 0:
            roots.append((e[0] + sum(comb[j]*e[j+1] for j in xrange(n-1)))/2)

    assert len(roots) == 120

    simple_roots = [
        (e[0]-e[1]-e[2]-e[3]-e[4]-e[5]-e[6]+e[7])/2,
        e[0]+e[1],
        e[0]-e[1],
        e[1]-e[2],
        e[2]-e[3],
        e[3]-e[4],
        e[4]-e[5],
        e[5]-e[6],
    ]
    highest_root = e[6] + e[7]

    assert all(alpha in roots for alpha in simple_roots + [highest_root])

    # restrict to get E6, E7, or E8
    simple_roots = simple_roots[:m]
    W = V.subspace(simple_roots)
    assert W.dimension() == m
    roots = [alpha for alpha in roots if alpha in W]

    if m == 7:
        # E7
        assert len(roots) == 63
        highest_root = e[6]-e[7]

    elif m == 6:
        # E6
        assert len(roots) == 36
        highest_root = (e[0]+e[1]+e[2]+e[3]+e[4]-e[5]-e[6]+e[7])/2

    assert highest_root in W

elif family == 'F':
    roots = [e[i] - e[j] for i in xrange(n) for j in xrange(i+1, n)] + \
        [e[i] + e[j] for i in xrange(n) for j in xrange(i+1, n)] + \
        [e[i] for i in xrange(n)] + \
        [(e[0] + sum(comb[j]*e[j+1] for j in xrange(n-1)))/2 for comb in itertools.product([1,-1], repeat=n-1)]

    assert len(roots) == 24

    simple_roots = [
        e[1]-e[2],
        e[2]-e[3],
        e[3],
        (e[0]-e[1]-e[2]-e[3])/2,
    ]
    highest_root = e[0]+e[1]

elif family == 'C':
    roots = [e[i] - e[j] for i in xrange(n) for j in xrange(i+1, n)] + \
        [e[i] + e[j] for i in xrange(n) for j in xrange(i+1, n)] + \
        [2*e[i] for i in xrange(n)]

    assert len(roots) == n**2

    simple_roots = [e[i]-e[i+1] for i in xrange(n-1)] + [2*e[n-1]]
    highest_root = 2*e[0]


assert all(alpha in roots for alpha in simple_roots + [highest_root])
assert len(simple_roots) == m

print "Simple roots:", simple_roots
print "Highest root:", highest_root
