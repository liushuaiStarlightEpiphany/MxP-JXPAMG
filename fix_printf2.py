f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

for i, line in enumerate(lines):
    if 'jx_printf([ClsFix]' in line:
        lines[i] = line.replace('[ClsFix]', '"[ClsFix]')
        lines[i] = lines[i].rstrip() + '\\n",\n'
        print(f'Fixed L{i+1}')
    if 'jx_printf( Proc %d: interp_type = 62' in line:
        lines[i] = '               jx_printf(" Proc %d: interp_type = 62 requires ai_measure_type = 1\\n", my_id);\n'
        print(f'Fixed L{i+1}')
    if lines[i].strip() == '",' and i > 0 and '\\n"' in lines[i-1]:
        lines[i] = ''

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
