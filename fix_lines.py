f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

for i, line in enumerate(lines):
    if 'jx_printf([EI-ideal]' in line:
        lines[i] = line.replace('[EI-ideal]', '"[EI-ideal]')
        lines[i] = lines[i].rstrip() + '\\n",\n'
        print(f'Fixed line {i+1}: {lines[i].strip()[:60]}')
    if 'ahi_level_num_ei_fpts = ahi_level_num_fpts;' in line and i > 1100:
        print(f'Check L{i+1}: {line.strip()}')

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
