f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

for i in range(len(lines)):
    # Fix any printf with missing quotes
    if 'jx_printf( Proc' in lines[i]:
        lines[i] = '               jx_printf(" Proc %d: interp_type = 61 requires ai_measure_type = 1\\n", my_id);\n'
        print(f'Fixed L{i+1}')
    # Remove any orphaned '",' lines
    if lines[i].strip() == '",':
        lines[i] = ''
        print(f'Removed orphan L{i+1}')

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
