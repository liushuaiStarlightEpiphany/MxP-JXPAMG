import os, math, random
random.seed(123)

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
        lam.sort(reverse=True); return lam
    lam = [a11, a22, a33]; lam.sort(reverse=True); return lam

def gen_diag(ratio):
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33); b = ratio * m
        a12 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a13 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a23 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        if a12!=0 and a13!=0 and a23!=0: break
    return [a11,a22,a33,a12,a13,a23]

def gen_ps1(ratio):
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33); b = ratio * m
        a12 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a13 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a23 = random.uniform(-1,1)*min(a22,a33)*10**random.uniform(-4,0)
        if a12!=0 and a13!=0 and a23!=0 and abs(a22)>0.5*(abs(a12)+abs(a23)) and abs(a33)>0.5*(abs(a13)+abs(a23)): break
    return [a11,a22,a33,a12,a13,a23]

def gen_ps2(ratio):
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33); b = ratio * m
        a12 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a23 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a13 = random.uniform(-1,1)*min(a11,a33)*10**random.uniform(-4,0)
        if a12!=0 and a13!=0 and a23!=0 and abs(a11)>0.5*(abs(a12)+abs(a13)) and abs(a33)>0.5*(abs(a13)+abs(a23)): break
    return [a11,a22,a33,a12,a13,a23]

def gen_ps3(ratio):
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33); b = ratio * m
        a13 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a23 = random.uniform(-1,1)*b*random.uniform(0.01,1)
        a12 = random.uniform(-1,1)*min(a11,a22)*10**random.uniform(-4,0)
        if a12!=0 and a13!=0 and a23!=0 and abs(a11)>0.5*(abs(a12)+abs(a13)) and abs(a22)>0.5*(abs(a12)+abs(a23)): break
    return [a11,a22,a33,a12,a13,a23]

def gen_general():
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33); ratio = 10**random.uniform(-4,0)
        a12 = random.uniform(-1,1)*m*ratio; a13 = random.uniform(-1,1)*m*ratio; a23 = random.uniform(-1,1)*m*ratio
        if a12!=0 and a13!=0 and a23!=0 and abs(a11)>0.5*(abs(a12)+abs(a13)) and abs(a22)>0.5*(abs(a12)+abs(a23)) and abs(a33)>0.5*(abs(a13)+abs(a23)): break
    return [a11,a22,a33,a12,a13,a23]

base = os.path.expanduser('~/l_s/eigenvalue-3/matrix-test/sweep')
ratios = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
names = ['1e-2', '1e-3', '1e-4', '1e-5', '1e-6']

for ratio, name in zip(ratios, names):
    d = os.path.join(base, name)
    os.makedirs(d, exist_ok=True)
    idx = 0
    for _ in range(20): write_matrix(d, idx, gen_diag(ratio), eigvals(gen_diag(ratio))); idx += 1
    for _ in range(20): write_matrix(d, idx, gen_ps1(ratio), eigvals(gen_ps1(ratio))); idx += 1
    for _ in range(20): write_matrix(d, idx, gen_ps2(ratio), eigvals(gen_ps2(ratio))); idx += 1
    for _ in range(20): write_matrix(d, idx, gen_ps3(ratio), eigvals(gen_ps3(ratio))); idx += 1
    for _ in range(20): write_matrix(d, idx, gen_general(), eigvals(gen_general())); idx += 1
    print('Done %s (%d)' % (name, idx))
print('All done')
