f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

for i in range(len(lines)):
    if 'interp_type = 61 requires' in lines[i]:
        lines[i] = '               jx_printf(" Proc %d: interp_type = 61 requires ai_measure_type = 1\\n", my_id);\n'
        # Remove the next line if it contains my_id
        if i+1 < len(lines) and 'my_id' in lines[i+1]:
            lines[i+1] = ''
        print(f'Fixed printf at L{i+1}')
        break

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
