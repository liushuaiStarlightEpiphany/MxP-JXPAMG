import os, math, random
random.seed(7)

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

def gen_ps1_dual():
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33)
        # row1 off-diags <= 1e-6 * m (trigger at eps=1e-6)
        a12 = random.uniform(-1,1) * 1e-6 * m * random.uniform(0.1,1.0)
        a13 = random.uniform(-1,1) * 1e-6 * m * random.uniform(0.1,1.0)
        # a23 between 1e-6 and 1e-4 (not trigger at 1e-6, but trigger at 1e-4)
        # so 1e-6*m < |a23| < 1e-4*m and satisfy diag dom
        scale = 10**random.uniform(-5.5, -4.1)
        a23 = random.uniform(-1,1) * m * scale
        if a12!=0 and a13!=0 and a23!=0 and abs(a23)/m > 1.2e-6:
            # Check diag dom
            if abs(a22)>0.5*(abs(a12)+abs(a23)) and abs(a33)>0.5*(abs(a13)+abs(a23)):
                break
    return [a11,a22,a33,a12,a13,a23]

def gen_ps2_dual():
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33)
        a12 = random.uniform(-1,1) * 1e-6 * m * random.uniform(0.1,1.0)
        a23 = random.uniform(-1,1) * 1e-6 * m * random.uniform(0.1,1.0)
        scale = 10**random.uniform(-5.5, -4.1)
        a13 = random.uniform(-1,1) * m * scale
        if a12!=0 and a13!=0 and a23!=0 and abs(a13)/m > 1.2e-6:
            if abs(a11)>0.5*(abs(a12)+abs(a13)) and abs(a33)>0.5*(abs(a13)+abs(a23)):
                break
    return [a11,a22,a33,a12,a13,a23]

def gen_ps3_dual():
    while True:
        a11 = 10**random.uniform(0, 12); a22 = 10**random.uniform(0, 12); a33 = 10**random.uniform(0, 12)
        m = min(a11, a22, a33)
        a13 = random.uniform(-1,1) * 1e-6 * m * random.uniform(0.1,1.0)
        a23 = random.uniform(-1,1) * 1e-6 * m * random.uniform(0.1,1.0)
        scale = 10**random.uniform(-5.5, -4.1)
        a12 = random.uniform(-1,1) * m * scale
        if a12!=0 and a13!=0 and a23!=0 and abs(a12)/m > 1.2e-6:
            if abs(a11)>0.5*(abs(a12)+abs(a13)) and abs(a22)>0.5*(abs(a12)+abs(a23)):
                break
    return [a11,a22,a33,a12,a13,a23]

base = os.path.expanduser('~/l_s/eigenvalue-3/matrix-test/dual')
os.makedirs(base, exist_ok=True)

# 100 ps1-like, 100 ps2-like, 100 ps3-like
for i in range(100):
    A = gen_ps1_dual(); write_matrix(base, i, A, eigvals(A))
for i in range(100, 200):
    A = gen_ps2_dual(); write_matrix(base, i, A, eigvals(A))
for i in range(200, 300):
    A = gen_ps3_dual(); write_matrix(base, i, A, eigvals(A))

print('Done 300 matrices')

# Verify
import os
t6 = 0; t4 = 0
for f in os.listdir(base):
    with open(os.path.join(base, f)) as fp:
        v = [float(fp.readline()) for _ in range(6)]
    a11,a22,a33,a12,a13,a23=v
    m = min(abs(a11),abs(a22),abs(a33))
    r12,r13,r23 = abs(a12)/m, abs(a13)/m, abs(a23)/m
    # At eps=1e-4: all three rows trigger?
    all_small_4 = (r12<1e-4 and r13<1e-4 and r23<1e-4)
    # At eps=1e-6: exactly one row triggers?
    r1_small_6 = (r12<1e-6 and r13<1e-6)
    r2_small_6 = (r12<1e-6 and r23<1e-6)
    r3_small_6 = (r13<1e-6 and r23<1e-6)
    one_row_6 = (r1_small_6 or r2_small_6 or r3_small_6) and not (r1_small_6 and r2_small_6)
    if one_row_6: t6 += 1
    if all_small_4: t4 += 1
print('Trigger at 1e-6 (one row): %d/300' % t6)
print('Trigger at 1e-4 (all rows): %d/300' % t4)
print('Both conditions met: %d/300' % (t6 if t6==t4==300 else 0))
