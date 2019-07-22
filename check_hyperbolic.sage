import argparse
import itertools
from operator import mul

import multiprocessing
# https://stackoverflow.com/questions/31941951/workaround-memory-leak-in-shared-object
def use_subprocess(func):
      def conn_func(conn, *args, **kwargs):
            conn.send(func(*args, **kwargs))
            conn.close()

      def new_function(*args, **kwargs):
            parent_conn, child_conn = multiprocessing.Pipe()
            p = multiprocessing.Process(target=conn_func, args=[child_conn]+list(args), kwargs=kwargs)
            p.start()
            result = parent_conn.recv()
            p.join()
            return result

      return new_function

###########################################################

parser = argparse.ArgumentParser(description='Check hyperbolic isometries of [1,w].')
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

e = V.basis()
assert all(e[i] == vector(matrix.identity(n)[i,:]) for i in xrange(n))


# construct root systems [Hum92, Section 2.10]
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
        [(e[0] + sum(comb[j]*e[j+1] for j in xrange(n-1)))/2
            for comb in itertools.product([1,-1], repeat=n-1)]

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

print "Simple roots: {}".format(simple_roots)
print "Highest root: {}".format(highest_root)

for alpha in roots:
    alpha.set_immutable()


def reflection(alpha, k=0):
    """
    Reflection w.r.t. the hyperplane {alpha*x = k}.
    """
    # find a point of the form t * alpha fixed by the reflection
    t = k / (alpha*alpha)
    p = t * alpha
    return Aff.translation(p) * Aff.reflection(alpha) * Aff.translation(p)**(-1)

# construct simple reflections [Hum92, Section 4.3]
S = [reflection(alpha, 0) for alpha in simple_roots] + [reflection(highest_root, 1)]


# Coxeter element
w = reduce(mul, S, Aff.one())

print "Coxeter element:"
print w

# find order of the linear part of w (so that w**order is a pure translation)
order = 1
while (w**order).matrix()[:n,:n] != matrix.identity(n):
    order += 1

# find minimal move vector of w (which is also the direction of the Coxeter axis)
axis_direction = (w**order)(zero_vector(n)) / order

# find a point on the Coxeter axis
A = w.matrix()[:n,:n]
b = vector(w.matrix()[:n,n])
point_on_axis = (A-A.parent().one()).solve_right(axis_direction - b)

assert w(point_on_axis) == point_on_axis + axis_direction

# make sure that it is in the interior of some chamber
if any(alpha*point_on_axis in ZZ for alpha in roots):
    # slide the point upwards
    point_on_axis += axis_direction / 4
assert all(alpha*point_on_axis not in ZZ for alpha in roots)

print "Coxeter axis: {} + theta * {}".format(point_on_axis, axis_direction)


# find horizontal roots
horizontal_roots = [alpha for alpha in roots if alpha*axis_direction == 0]
print "Number of pairs of horizontal roots: {}".format(len(horizontal_roots))

# split horizontal roots into irreducible components
horizontal_components = DisjointSet(horizontal_roots)

for alpha, beta in itertools.combinations(horizontal_roots, 2):
    if alpha*beta != 0:
        horizontal_components.union(alpha, beta)

horizontal_components = list(sorted(horizontal_components, key=len))
sizes = [V.subspace(component).dimension() for component in horizontal_components]
print "Horizontal components: {}".format(', '.join(["A{}".format(k) for k in sizes]))



# find horizontal elements in [1,w]
horizontal_elements_by_component = [[] for i in xrange(len(horizontal_components))]

for i in xrange(len(horizontal_components)):
    size = sizes[i]
    component = horizontal_components[i]
    W = V.subspace(component)   # subspace of V generated by the roots in this component
    basis = W.basis()

    # find elements that fix the origin
    elements = []

    for comb in itertools.permutations(component, size):
        wh = reduce(mul, [reflection(alpha) for alpha in comb], Aff.one())

        if all(w(v) - w(zero_vector(n)) == wh(v) for v in basis):
            # the linear part of wh agrees with the linear part of w on this irreducible component,
            # so wh is the (unique) maximal element of this component that fixes the origin

            for j in xrange(size+1):
                elements.append(reduce(mul, [reflection(alpha) for alpha in comb[:j]], Aff.one()))

    # remove duplicates
    elements = [element for element, group in itertools.groupby(sorted(elements))]

    # find all elements by conjugation by powers of w
    for u in elements:
        for j in xrange(size+1):
            horizontal_elements_by_component[i].append((w**j) * u * (w**(-j)))

    # remove duplicates
    horizontal_elements_by_component[i] = [
        element for element, group in itertools.groupby(sorted(horizontal_elements_by_component[i]))
    ]

    # a horizontal component should form half of a noncrossing partition lattice of type B_{size+1}
    # and the cardinality of a noncrossing partition lattice is given in [Arm09, Figure 2.8]
    assert len(horizontal_elements_by_component[i]) == binomial(2*(size+1), size+1) / 2

horizontal_elements = [reduce(mul, comb, Aff.one()) for comb in itertools.product(*horizontal_elements_by_component)]


def length(u):
    """
    Reflection length of an elliptic isometry.
    """
    A = u.matrix()[:n,:n]   # linear part of u
    return n - (A - A.parent().one()).kernel().dimension()

# sort horizontal elements by reflection length
horizontal_elements = sorted(horizontal_elements, key=lambda u: length(u))

print "Horizontal elements by length: {}".format([
        len(list(group)) for l, group in itertools.groupby(horizontal_elements, key=lambda u: length(u))
    ])  # in the case E8, this is equal to the bottom row of [MS17, Figure 13]



# check hyperbolic elements
print "\nCheck hyperbolic elements..."

orthogonal = True   # are positive (resp. negative) walls always pairwise orthogonal?

@use_subprocess
def get_walls(n, p, hyperplanes):
    """
    Given a list of hyperplanes {alpha*x = k}, given as pairs (alpha, k), and a
    point p in the complement of every hyperplane, find the walls of the chamber
    containing p.
    """
    chamber = Polyhedron(ieqs=[
        ([-k] + list(alpha) if alpha*p - k > 0 else [k] + list(-alpha)) \
        for alpha, k in hyperplanes])
    walls = []

    for face in chamber.faces(n-1):
        [equation] = face.as_polyhedron().equations()
        k = -equation[0]
        alpha = vector(equation[1:])

        walls.append((alpha, k))

    return walls

for i, h in enumerate(horizontal_elements):
    u = w*(h**(-1)) # left complement of the horizontal element h
    assert u*h == w
    l = m+1-length(h)   # reflection length of u
    print "[{}/{}]".format(i+1, len(horizontal_elements))

    if args.verbose:
        print "l(u) = {}".format(l)
        print u
        print

    # find W = Span(Mov(u))
    W = V.subspace([V(u(p)-p) for p in e + [zero_vector(n)]])
    assert W.dimension() == l-1

    # find order of the linear part of u
    u_order = 1
    A = u.matrix()[:n,:n]
    while A**u_order != matrix.identity(n):
        u_order += 1
    assert (u**u_order).matrix()[:n,:n] == matrix.identity(n)

    # find roots below u (necessary condition)
    available_roots = [alpha for alpha in roots if alpha in W]
    available_horizontal_roots = [alpha for alpha in available_roots if alpha*axis_direction == 0]
    available_vertical_roots = [alpha for alpha in available_roots if alpha*axis_direction != 0]

    available_horizontal_hyperplanes = []

    for alpha in available_horizontal_roots:
        # which horizontal hyperplanes {alpha*x = k} are really below u?
        # r <= u if and only if r*h is a horizontal element of [1,w] and l(r*h)=l(h)+1
        ok = False
        for k in [floor(alpha*point_on_axis), floor(alpha*point_on_axis)+1]:
            r = reflection(alpha, k)
            if r*h in horizontal_elements and length(r*h) == length(h)+1:
                available_horizontal_hyperplanes.append((alpha, k))
                ok = True

        assert ok   # there should be at least one horizontal reflection below u

    for j in xrange(2*order):
        a = point_on_axis + j * axis_direction/2

        # p is not on any hyperplane
        assert all(alpha*a not in ZZ for alpha in roots)

        # find possible walls of p
        available_hyperplanes = list(itertools.chain.from_iterable(
            [(alpha, floor(alpha*a)), (alpha, floor(alpha*a)+1)] for alpha in available_vertical_roots)) \
            + available_horizontal_hyperplanes

        walls = get_walls(n, a, available_hyperplanes)
        assert len(walls) == l

        positive_walls = sorted(
            [(alpha, k) for (alpha, k) in walls if alpha*axis_direction/(k-alpha*a) > 0],
            key=lambda (alpha, k): (k-alpha*a)/(alpha*axis_direction))
        negative_walls = sorted(
            [(alpha, k) for (alpha, k) in walls if alpha*axis_direction/(k-alpha*a) < 0],
            key=lambda (alpha, k): (k-alpha*a)/(alpha*axis_direction))
        horizontal_walls = [(alpha, k) for (alpha, k) in walls if alpha*axis_direction == 0]

        if args.verbose and orthogonal:
            # check if positive walls are orthogonal
            if any(alpha*beta != 0 for (alpha, k), (beta, h) in itertools.combinations(positive_walls, 2)):
                print "The positive/negative walls are not pairwise orthogonal!"
                orthogonal = False

            # check if negative walls are orthogonal
            if any(alpha*beta != 0 for (alpha, k), (beta, h) in itertools.combinations(negative_walls, 2)):
                print "The positive/negative walls are not pairwise orthogonal!"
                orthogonal = False

        # check that u can be obtained as the product of the walls
        assert any(reduce(mul, [reflection(alpha, k) for alpha, k in positive_walls + list(horizontal) + negative_walls], Aff.one()) == u \
            for horizontal in itertools.permutations(horizontal_walls))

print "Case {} checked successfully".format(args.case)
