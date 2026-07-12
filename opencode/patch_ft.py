#!/usr/bin/env python3
"""Modify ft_test.c and ft_compare.c to use robust for Case 5 in a2()."""

import sys

def patch_file(fpath):
    with open(fpath, 'r') as f:
        content = f.read()
    
    # 1. Replace Cardano code in a2() Case 5 with robust call
    # The Cardano code starts with "double t3 = (m[0] + m[1] + m[2]) / 3;"
    old_cardano_block = """    double t3 = (m[0] + m[1] + m[2]) / 3;
    double p = m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
             - m[3]*m[3] - m[4]*m[4] - m[5]*m[5] - 3*t3*t3;
    double q = m[0]*m[1]*m[2] + 2*m[3]*m[5]*m[4]
             - m[0]*m[5]*m[5] - m[1]*m[4]*m[4] - m[2]*m[3]*m[3]
             + 2*t3*t3*t3 - t3*(m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
             - m[3]*m[3] - m[4]*m[4] - m[5]*m[5]);
    double dp = p/3, dq = -q/2, dis = dq*dq + dp*dp*dp;
    if (dis >= 0) {
        double s = sqrt(dis), u = cbrt(dq + s), v = cbrt(dq - s);
        l[0] = u+v+t3; l[1] = -(u+v)/2+t3; l[2] = l[1];
        l[1] += sqrt(3)/2*(u-v); l[2] -= sqrt(3)/2*(u-v);
    } else {
        double phi = acos(dq/sqrt(-dp*dp*dp))/3, r = 2*sqrt(-dp);
        l[0] = r*cos(phi)+t3; l[1] = r*cos(phi+2*M_PI/3)+t3; l[2] = r*cos(phi+4*M_PI/3)+t3;
    }
    qsort(l, 3, sizeof(double), cd);"""

    new_robust_block = """    double rr[6];
    rr[0] = m[0]; rr[1] = m[3]; rr[2] = m[4];
    rr[3] = m[1]; rr[4] = m[5]; rr[5] = m[2];
    robust(rr, l);
    qsort(l, 3, sizeof(double), cd);"""

    if old_cardano_block in content:
        content = content.replace(old_cardano_block, new_robust_block)
        print(f"  Replaced Cardano block in a2()")
    else:
        # Try with PI instead of M_PI
        old_cardano_block2 = old_cardano_block.replace("M_PI", "PI")
        if old_cardano_block2 in content:
            content = content.replace(old_cardano_block2, new_robust_block)
            print(f"  Replaced Cardano block in a2() (PI variant)")
        else:
            print(f"  WARNING: Could not find Cardano block in a2()!")
            print(f"  File snippet:")
            # Find where double t3 is
            idx = content.find("double t3 = ")
            if idx >= 0:
                print(content[idx:idx+600])
            return False

    # 2. Add cardano_a2 function (insert before a2() or after it)
    cardano_func = """
static void cardano_a2(const double m[6], double l[3])
{
    double t3 = (m[0] + m[1] + m[2]) / 3;
    double p = m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
             - m[3]*m[3] - m[4]*m[4] - m[5]*m[5] - 3*t3*t3;
    double q = m[0]*m[1]*m[2] + 2*m[3]*m[5]*m[4]
             - m[0]*m[5]*m[5] - m[1]*m[4]*m[4] - m[2]*m[3]*m[3]
             + 2*t3*t3*t3 - t3*(m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
             - m[3]*m[3] - m[4]*m[4] - m[5]*m[5]);
    double dp = p/3, dq = -q/2, dis = dq*dq + dp*dp*dp;
    if (dis >= 0) {
        double s = sqrt(dis), u = cbrt(dq + s), v = cbrt(dq - s);
        l[0] = u+v+t3; l[1] = -(u+v)/2+t3; l[2] = l[1];
        l[1] += sqrt(3)/2*(u-v); l[2] -= sqrt(3)/2*(u-v);
    } else {
        double phi = acos(dq/sqrt(-dp*dp*dp))/3, r = 2*sqrt(-dp);
        l[0] = r*cos(phi)+t3; l[1] = r*cos(phi+2*M_PI/3)+t3; l[2] = r*cos(phi+4*M_PI/3)+t3;
    }
    qsort(l, 3, sizeof(double), cd);
}
"""
    # Insert cardano_a2 after the a2 function
    # Find end of a2 function - after qsort(l, 3, sizeof(double), cd);} followed by static
    end_of_a2 = content.find("qsort(l, 3, sizeof(double), cd);\n}")
    if end_of_a2 >= 0:
        end_of_a2 = content.find("\n}", end_of_a2) + 2
        content = content[:end_of_a2] + cardano_func + content[end_of_a2:]
        print(f"  Added cardano_a2() after a2()")
    else:
        print(f"  WARNING: Could not find end of a2()!")
        return False

    # 3. Replace a2(m, l, eps_zero_pattern) with cardano_a2(m, l)
    # Look for a2(m, l, 0) - standalone calls
    import re
    # Count occurrences
    calls_a2_zero = [(m.start(), m.group()) for m in re.finditer(r'a2\([^)]*,\s*0\s*\)', content)]
    calls_printed = set()
    
    for start, match in calls_a2_zero:
        if match not in calls_printed:
            calls_printed.add(match)
            print(f"  Found: {match}")
    
    # Replace a2(m, l, 0) with cardano_a2(m, l)
    content = re.sub(r'a2\(([^,]+),\s*([^,]+),\s*0\s*\)', r'cardano_a2(\1, \2)', content)
    print(f"  Replaced a2(*, *, 0) with cardano_a2(*, *)")

    # 4. Check if M_PI is used (need math.h) - already included

    with open(fpath, 'w') as f:
        f.write(content)
    
    print(f"  Written to {fpath}")
    return True

if __name__ == '__main__':
    for f in sys.argv[1:]:
        print(f"Patching {f}...")
        if patch_file(f):
            print(f"  OK")
        else:
            print(f"  FAILED")
            sys.exit(1)
