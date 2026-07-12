with open("/vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/src/test_bsr_cprbicgstab_hybrid.c", "r") as f:
    content = f.read()

# Find the solve function and fix the CPR setup issue
old = '''    // 使用低精度矩阵设置CPR
    double setup_start = MPI_Wtime();
    JXF_Int setup_result = JXF_CPRSetup(cpr, A_low);
    double setup_end = MPI_Wtime();
    
    if (setup_result != JXF_SUCCESS) {
        if (myid == 0) printf("Error setting up low precision CPR preconditioner\\n");
        jxf_ParBSRMatrixDestroy(A_low);
        JXF_CPRDestroy(&cpr);
        return -1;
    }
    
    if (myid == 0) {
        printf("Low precision CPR setup completed in %.6f seconds\\n", setup_end - setup_start);
    }
    
    // 3. 创建高精度BiCGSTAB求解器
    // BSR BiCGSTAB
    JX_Solver bicgstab_solver;
    jx_BiCGSTABFunctions *fg = jx_BiCGSTABFunctionsCreate(
        jx_ParBSRKrylovCreateVector, jx_ParBSRKrylovDestroyVector,
        jx_ParBSRKrylovMatvecCreate, jx_ParBSRKrylovMatvec, jx_ParBSRKrylovMatvecDestroy,
        jx_ParBSRKrylovInnerProd, jx_ParBSRKrylovCopyVector,
        jx_ParBSRKrylovClearVector, jx_ParBSRKrylovScaleVector, jx_ParBSRKrylovAxpy,
        jx_ParBSRKrylovCommInfo,
        (JX_Int (*)(void*, void*))JXF_CPRSetup,
        (JX_Int (*)(void*, void*, void*, void*))JXF_MxP_CPRPrecond1);
    bicgstab_solver = (JX_Solver)jx_BiCGSTABCreate(fg);'''

new = '''    // 创建一个包装函数，使得 BiCGSTAB 的 setup 阶段使用低精度矩阵
    // 这里的 cpr_ptr 指向 cpr，matrix_ptr 是 BiCGSTAB 传入的 A_high（忽略）
    // 我们使用前面已经创建好的 A_low
    double setup_start = MPI_Wtime();
    JXF_Int setup_result = JXF_CPRSetup(cpr, A_low);
    double setup_end = MPI_Wtime();
    
    if (setup_result != JXF_SUCCESS) {
        if (myid == 0) printf("Error setting up low precision CPR preconditioner\\n");
        jxf_ParBSRMatrixDestroy(A_low);
        JXF_CPRDestroy(&cpr);
        return -1;
    }
    
    if (myid == 0) {
        printf("Low precision CPR setup completed in %.6f seconds\\n", setup_end - setup_start);
    }
    
    // 3. 创建高精度BiCGSTAB求解器
    // BSR BiCGSTAB
    JX_Solver bicgstab_solver;
    jx_BiCGSTABFunctions *fg = jx_BiCGSTABFunctionsCreate(
        jx_ParBSRKrylovCreateVector, jx_ParBSRKrylovDestroyVector,
        jx_ParBSRKrylovMatvecCreate, jx_ParBSRKrylovMatvec, jx_ParBSRKrylovMatvecDestroy,
        jx_ParBSRKrylovInnerProd, jx_ParBSRKrylovCopyVector,
        jx_ParBSRKrylovClearVector, jx_ParBSRKrylovScaleVector, jx_ParBSRKrylovAxpy,
        jx_ParBSRKrylovCommInfo,
        (JX_Int (*)(void*, void*))JXF_CPRSkipSetup,
        (JX_Int (*)(void*, void*, void*, void*))JXF_MxP_CPRPrecond1);
    bicgstab_solver = (JX_Solver)jx_BiCGSTABCreate(fg);'''

if old in content:
    content = content.replace(old, new)
    with open("/vol8/home/xtu_pcy/l_s/jxpamg_all/JXFPANG-bsr/src/test_bsr_cprbicgstab_hybrid.c", "w") as f:
        f.write(content)
    print("FIXED: CPR setup wrapper")
else:
    print("OLD TEXT NOT FOUND")
    idx = content.find("使用低精度矩阵设置CPR")
    if idx >= 0:
        print(f"Found at {idx}")
        print(repr(content[idx:idx+500]))
