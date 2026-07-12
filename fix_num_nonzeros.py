f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_interp_clsfix.c', 'r')
lines = f.readlines()
f.close()

for i in range(len(lines)):
    # After diag realloc, add num_nonzeros update
    if 'for (JX_Int k = nz - 1; k >= cde; k--)' in lines[i]:
        # The next few lines do: shift, then update row pointers
        # After the row pointer update loop, add:
        for j in range(i, i+20):
            if 'for (JX_Int k = i+1; k <= nr; k++) cdi[k] += ediff;' in lines[j]:
                lines.insert(j+1, '               jx_CSRMatrixNumNonzeros(cd) = nz + ediff;\n')
                print(f'Added diag num_nonzeros update after line {j+1}')
                break
        break

# Similarly for offd
for i in range(len(lines)):
    if 'for (JX_Int k = nz - 1; k >= coe; k--)' in lines[i]:
        for j in range(i, i+20):
            if 'for (JX_Int k = i+1; k <= jx_CSRMatrixNumRows(co); k++) coi[k] += odiff;' in lines[j]:
                lines.insert(j+1, '               jx_CSRMatrixNumNonzeros(co) = nz + odiff;\n')
                print(f'Added offd num_nonzeros update after line {j+1}')
                break
        break

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/csrc/amg/par_interp_clsfix.c', 'w')
f.writelines(lines)
f.close()
print('Done')
