f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'r')
lines = f.readlines()
f.close()

# Remove the orphaned ',' line
for i in range(len(lines)):
    if lines[i].strip() == ',' and i > 0 and '\\n",' in lines[i-1]:
        lines.pop(i)
        print(f'Removed orphaned comma at line {i+1}')
        break

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_amg_setup.c', 'w')
f.writelines(lines)
f.close()
print('Done')
