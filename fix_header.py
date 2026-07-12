f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/include/jx_pamg.h', 'r')
lines = f.readlines()
f.close()

# Find the end of PAMGBuildExtPIInterp_mix declaration
for i, line in enumerate(lines):
    if 'PAMGBuildExtPIInterp_mix' in line and 'JX_Int' in lines[i-1]:
        # Find the closing semicolon
        for j in range(i, i+15):
            if ');' in lines[j]:
                # Add new declaration after this
                new_lines = [
                    lines[j].replace(');', ');\n'),
                    'JX_Int\n',
                    'jx_PAMGBuildClassicalFix(jx_ParCSRMatrix *par_A,\n',
                    '                       JX_Int *CF_marker,\n',
                    '                       jx_ParCSRMatrix *par_S,\n',
                    '                       JX_Real *eta,\n',
                    '                       JX_Real eta_threshold,\n',
                    '                       JX_Int *num_cpts_global,\n',
                    '                       JX_Int num_functions,\n',
                    '                       JX_Int *dof_func,\n',
                    '                       JX_Int debug_flag,\n',
                    '                       JX_Real trunc_factor,\n',
                    '                       JX_Int max_elmts,\n',
                    '                       JX_Int *col_offd_S_to_A,\n',
                    '                       jx_ParCSRMatrix **P_ptr,\n',
                    '                       JX_Int *num_fpts_global_ptr,\n',
                    '                       JX_Int *num_ei_fpts_global_ptr);\n',
                ]
                lines[j] = ''.join(new_lines)
                print(f'Added declaration after line {j+1}')
                break
        break

f = open('/vol8/home/xtu_pcy/l_s/jxpamg-basic-ipt/include/jx_pamg.h', 'w')
f.writelines(lines)
f.close()
print('Done')
