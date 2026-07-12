import os, math, random

random.seed(123)

def eigvals(A):
    a11, a22, a33, a12, a13, a23 = A
    tr = (a11 + a22 + a33) / 3.0
    A0a, A0b, A0c = a11-tr, a22-tr, a33-tr
    J2 = a12**2 + a13**2 + a23**2 - A0a*A0b - A0a*A0c - A0b*A0c
    J3 = A0a*A0b*A0c + a12*a23*a13 + a13*a12*a23 - a13*A0b*a13 - a23*a23*A0a - a12*a12*A0c
    if J2 > 0:
        t = max(-1.0, min(1.0, J3 / (2.0 * math.sqrt(J2**3 / 27.0))))
        alpha = math.acos(t) / 3.0
        s = 2.0 * math.sqrt(J2 / 3.0)
        lam = [s*math.cos(alpha)+tr, s*math.cos(alpha+2*math.pi/3)+tr, s*math.cos(alpha+4*math.pi/3)+tr]
        lam.sort(reverse=True)
        return lam
    lam = [a11, a22, a33]
    lam.sort(reverse=True)
    return lam

base = os.path.expanduser('~/l_s/eigenvalue-3/matrix-test/sweep')
ratios = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
names = ['1e-2', '1e-3', '1e-4', '1e-5', '1e-6']

for name in names:
    os.makedirs(os.path.join(base, name), exist_ok=True)

for ratio, name in zip(ratios, names):
    for i in range(100):
        a11 = 10**random.uniform(0, 12)
        a22 = 10**random.uniform(0, 12)
        a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33)
        bound = ratio * m
        r12 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
        r13 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
        r23 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
        A = [a11, a22, a33, r12, r13, r23]
        lam = eigvals(A)
        fname = os.path.join(base, name, 'matrix_%05d.txt' % i)
        with open(fname, 'w') as f:
            for v in A: f.write('%.15e\n' % v)
            for v in lam: f.write('%.15e\n' % v)
    print('Done %s' % name)
print('All done')
