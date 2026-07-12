import os, math, random

random.seed(42)

def write_matrix(outdir, idx, A, lam):
    fname = os.path.join(outdir, 'matrix_%05d.txt' % idx)
    with open(fname, 'w') as f:
        for v in A: f.write('%.15e\n' % v)
        for v in lam: f.write('%.15e\n' % v)

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

def gen_row_12():
    a11 = 10**random.uniform(0, 12)
    a22 = 10**random.uniform(0, 12)
    a33 = 10**random.uniform(0, 12)
    m = min(a11, a22, a33)
    # a12, a13 are small (controlled)
    bound = 5e-5 * m
    a12 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
    a13 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
    # a23 just satisfies diag dom: |a22| > 0.5*(|a12|+|a23|) and |a33| > 0.5*(|a13|+|a23|)
    # So |a23| < 2*min(|a22|,|a33|) - (other terms)
    while True:
        a23 = random.uniform(-1, 1) * min(a22, a33) * 10**random.uniform(-4, 0)
        if (abs(a22) > 0.5*(abs(a12)+abs(a23)) and 
            abs(a33) > 0.5*(abs(a13)+abs(a23)) and
            a23 != 0):
            break
    return [a11, a22, a33, a12, a13, a23]

def gen_row_13():
    a11 = 10**random.uniform(0, 12)
    a22 = 10**random.uniform(0, 12)
    a33 = 10**random.uniform(0, 12)
    m = min(a11, a22, a33)
    bound = 5e-5 * m
    a12 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
    a23 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
    while True:
        a13 = random.uniform(-1, 1) * min(a11, a33) * 10**random.uniform(-4, 0)
        if (abs(a11) > 0.5*(abs(a12)+abs(a13)) and
            abs(a33) > 0.5*(abs(a13)+abs(a23)) and
            a13 != 0):
            break
    return [a11, a22, a33, a12, a13, a23]

def gen_row_23():
    a11 = 10**random.uniform(0, 12)
    a22 = 10**random.uniform(0, 12)
    a33 = 10**random.uniform(0, 12)
    m = min(a11, a22, a33)
    bound = 5e-5 * m
    a13 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
    a23 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
    while True:
        a12 = random.uniform(-1, 1) * min(a11, a22) * 10**random.uniform(-4, 0)
        if (abs(a11) > 0.5*(abs(a12)+abs(a13)) and
            abs(a22) > 0.5*(abs(a12)+abs(a23)) and
            a12 != 0):
            break
    return [a11, a22, a33, a12, a13, a23]

def gen_general():
    while True:
        a11 = 10**random.uniform(0, 12)
        a22 = 10**random.uniform(0, 12)
        a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33)
        ratio = 10**random.uniform(-4, 0)
        a12 = random.uniform(-1, 1) * m * ratio
        a13 = random.uniform(-1, 1) * m * ratio
        a23 = random.uniform(-1, 1) * m * ratio
        if (abs(a11) > 0.5*(abs(a12)+abs(a13)) and
            abs(a22) > 0.5*(abs(a12)+abs(a23)) and
            abs(a33) > 0.5*(abs(a13)+abs(a23)) and
            a12 != 0 and a13 != 0 and a23 != 0):
            break
    return [a11, a22, a33, a12, a13, a23]

base = os.path.expanduser('~/l_s/eigenvalue-3/matrix-test')

# diag: all three off-diagonals <= 5e-5 * min(all diags)
for i in range(500):
    while True:
        a11 = 10**random.uniform(0, 12)
        a22 = 10**random.uniform(0, 12)
        a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33)
        bound = 5e-5 * m
        a12 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
        a13 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
        a23 = random.uniform(-1, 1) * bound * random.uniform(0.01, 1.0)
        if a12 != 0 and a13 != 0 and a23 != 0:
            break
    A = [a11, a22, a33, a12, a13, a23]
    write_matrix(os.path.join(base, 'diag'), i, A, eigvals(A))

# ps1: row1 off-diags <= 5e-5*min(all); others satisfy diag dom
for i in range(500):
    A = gen_row_12()
    write_matrix(os.path.join(base, 'ps1'), i, A, eigvals(A))

# ps2: row2 off-diags <= 5e-5*min(all); others satisfy diag dom
for i in range(500):
    A = gen_row_13()
    write_matrix(os.path.join(base, 'ps2'), i, A, eigvals(A))

# ps3: row3 off-diags <= 5e-5*min(all); others satisfy diag dom
for i in range(500):
    A = gen_row_23()
    write_matrix(os.path.join(base, 'ps3'), i, A, eigvals(A))

# general: all satisfy diag dom
for i in range(500):
    A = gen_general()
    write_matrix(os.path.join(base, 'general'), i, A, eigvals(A))

print('Done')
