f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

for i in range(len(lines)):
    # Remove the orphaned ',' on a line by itself before the printf args
    if i+1 < len(lines) and lines[i].strip() == ',' and '>= 0) ahi_level_ratio_F' in lines[i-1] if i > 0 else False:
        # Check if this is the orphaned comma
        if i+2 < len(lines) and 'level, fine_size,' in lines[i+2]:
            lines[i] = ''
            print(f'Removed orphan comma L{i+1}')
    # Fix duplicate my_id
    if ', my_id);\n, my_id);' in lines[i]:
        print(f'Fixing duplicate my_id L{i+1}')
        lines[i] = lines[i].replace(', my_id);\n, my_id);', ', my_id);\n')

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
