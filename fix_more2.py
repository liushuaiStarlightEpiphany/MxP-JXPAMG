f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

# Check around the suspected problematic area
for i in range(1137, 1145):
    print(f'L{i+1}: {repr(lines[i])}')
