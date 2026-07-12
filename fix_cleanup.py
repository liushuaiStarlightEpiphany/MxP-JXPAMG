with open("/vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/src/test_bsr_cprbicgstab_hybrid.c", "r") as f:
    content = f.read()

# Fix the cleanup: the A_low shares row_starts/col_starts with A_high,
# so after destroying A_low, the pointers in A_parbsr might be invalid.
# We need to destroy A_parbsr before A_low, or set shared pointers to NULL.

old_cleanup = '''    // 8. 清理
    jx_BiCGSTABDestroy(bicgstab_solver);
    JXF_CPRDestroy(&cpr);
    jxf_ParBSRMatrixDestroy(A_low);
    
    return 0;'''

new_cleanup = '''    // 8. 清理
    jx_BiCGSTABDestroy(bicgstab_solver);
    JXF_CPRDestroy(&cpr);
    /* A_low shares row_starts/col_starts with A_high; detach them before destroy */
    if (A_low) {
        A_low->row_starts = NULL;
        A_low->col_starts = NULL;
        A_low->diag = NULL;
        A_low->offd = NULL;
        A_low->comm = NULL;
    }
    jxf_ParBSRMatrixDestroy(A_low);
    
    return 0;'''

if old_cleanup in content:
    content = content.replace(old_cleanup, new_cleanup)
    with open("/vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/src/test_bsr_cprbicgstab_hybrid.c", "w") as f:
        f.write(content)
    print("FIXED: cleanup")
else:
    print("NOT FOUND")
    idx = content.find("8. 清理")
    if idx >= 0:
        print(repr(content[idx:idx+200]))
