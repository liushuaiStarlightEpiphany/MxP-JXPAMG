with open("/vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/src/test_bsr_cprbicgstab_hybrid.c", "r") as f:
    content = f.read()

# Add a no-op CPR setup wrapper function before solve_with_mixed_precision_cpr_gmres
old_func_start = "int solve_with_mixed_precision_cpr_gmres("
if old_func_start in content:
    wrapper = '''/* CPR setup wrapper: skip re-setup since CPR is already configured with A_low */
static JX_Int mixed_cpr_setup_wrapper(void *cpr, void *matrix)
{
    (void)matrix; /* ignore - A_low was already used for setup */
    return JXF_SUCCESS;
}

int solve_with_mixed_precision_cpr_gmres(
'''
    content = content.replace(old_func_start, wrapper)
    with open("/vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/src/test_bsr_cprbicgstab_hybrid.c", "w") as f:
        f.write(content)
    print("FIXED: added wrapper")
else:
    print("NOT FOUND")
