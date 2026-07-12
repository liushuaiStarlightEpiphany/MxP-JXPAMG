f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

for i in range(len(lines)):
    if 'jx_printf("[EI-ideal]' in lines[i]:
        # Check if next line is just '",'
        if i+1 < len(lines) and lines[i+1].strip() == '",':
            lines[i+1] = ''
            print(f'Removed orphaned line {i+2}')
    if 'jx_printf( Proc %d:' in lines[i]:
        lines[i] = lines[i].replace('jx_printf( Proc %d:', 'jx_printf(" Proc %d:')
        # Find the orphaned end
        if i+1 < len(lines) and lines[i+1].strip() == '",':
            lines[i+1] = ''
        print(f'Fixed printf at line {i+1}')

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
