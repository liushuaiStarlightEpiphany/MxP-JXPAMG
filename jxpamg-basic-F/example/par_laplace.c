//========================================================================//
//  JXFPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_laplace.c -- Laplace problem
 *  Date: 2012/09/19
 *
 *  Created by Yue Xiaoqiang
 */

#include "jxf_util.h"
#include "jxf_mv.h"

/*!
 * \fn JXF_Int jxf_map
 */
JXF_Int
jxf_map( JXF_Int ix,
        JXF_Int iy,
        JXF_Int iz,
        JXF_Int p,
        JXF_Int q,
        JXF_Int r,
        JXF_Int P,
        JXF_Int Q,
        JXF_Int R,
        JXF_Int *nx_part,
        JXF_Int *ny_part,
        JXF_Int *nz_part,
        JXF_Int *global_part )
{
    JXF_Int nx_local;
    JXF_Int ny_local;
    JXF_Int ix_local;
    JXF_Int iy_local;
    JXF_Int iz_local;
    JXF_Int global_index;
    JXF_Int proc_num;
    
    proc_num = r * P * Q + q * P + p;
    nx_local = nx_part[p+1] - nx_part[p];
    ny_local = ny_part[q+1] - ny_part[q];
    ix_local = ix - nx_part[p];
    iy_local = iy - ny_part[q];
    iz_local = iz - nz_part[r];
    global_index = global_part[proc_num] + (iz_local * ny_local + iy_local) * nx_local + ix_local;
    
    return global_index;
}

/*!
 * \fn JXF_Int jxf_map2
 */
JXF_Int
jxf_map2( JXF_Int ix,
         JXF_Int iy,
         JXF_Int p,
         JXF_Int q,
         JXF_Int P,
         JXF_Int Q,
         JXF_Int *nx_part,
         JXF_Int *ny_part,
         JXF_Int *global_part )
{
    JXF_Int nx_local;
    JXF_Int ix_local;
    JXF_Int iy_local;
    JXF_Int global_index;
    JXF_Int proc_num;
    
    proc_num = q * P + p;
    nx_local = nx_part[p+1] - nx_part[p];
    ix_local = ix - nx_part[p];
    iy_local = iy - ny_part[q];
    global_index = global_part[proc_num] + iy_local * nx_local + ix_local;
    
    return global_index;
}

/*!
 * \fn JXF_Int jxf_map3
 */
JXF_Int
jxf_map3( JXF_Int ix,
         JXF_Int iy,
         JXF_Int iz,
         JXF_Int p,
         JXF_Int q,
         JXF_Int r,
         JXF_Int P,
         JXF_Int Q,
         JXF_Int R,
         JXF_Int *nx_part,
         JXF_Int *ny_part,
         JXF_Int *nz_part,
         JXF_Int *global_part )
{
    JXF_Int nx_local;
    JXF_Int ix_local;
    JXF_Int iy_local;
    JXF_Int iz_local;
    JXF_Int nxy;
    JXF_Int global_index;
    JXF_Int proc_num;
    
    proc_num = r * P * Q + q * P + p;
    nx_local = nx_part[p+1] - nx_part[p];
    nxy = nx_local * (ny_part[q+1] - ny_part[q]);
    ix_local = ix - nx_part[p];
    iy_local = iy - ny_part[q];
    iz_local = iz - nz_part[r];
    global_index = global_part[proc_num] + iz_local * nxy + iy_local * nx_local + ix_local;
    
    return global_index;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_GenerateParLaplacian
 */
jxf_ParCSRMatrix *
jxf_GenerateParLaplacian( MPI_Comm comm,
                         JXF_Int nx,
                         JXF_Int ny,
                         JXF_Int nz,
                         JXF_Int P,
                         JXF_Int Q,
                         JXF_Int R,
                         JXF_Int p,
                         JXF_Int q,
                         JXF_Int r,
                         JXF_Real *value )
{
    jxf_ParCSRMatrix *A = NULL;
    
    jxf_CSRMatrix *diag = NULL;
    jxf_CSRMatrix *offd = NULL;
    
    JXF_Int *diag_i = NULL;
    JXF_Int *diag_j = NULL;
    JXF_Real *diag_data = NULL;
    
    JXF_Int *offd_i = NULL;
    JXF_Int *offd_j = NULL;
    JXF_Real *offd_data = NULL;
    
    JXF_Int *global_part = NULL;
    JXF_Int *col_map_offd = NULL;
    
    JXF_Int *nx_part = NULL;
    JXF_Int *ny_part = NULL;
    JXF_Int *nz_part = NULL;
    
    JXF_Int ix, iy, iz;
    JXF_Int cnt, o_cnt;
    JXF_Int local_num_rows;
    JXF_Int row_index;
    JXF_Int i, j;
    
    JXF_Int nx_local, ny_local, nz_local;
    JXF_Int nx_size, ny_size, nz_size;
    JXF_Int num_cols_offd;
    JXF_Int grid_size;
    
    JXF_Int num_procs, my_id;
    JXF_Int P_busy, Q_busy, R_busy;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny * nz;
    
    jxf_GeneratePartitioning(nx, P, &nx_part);
    jxf_GeneratePartitioning(ny, Q, &ny_part);
    jxf_GeneratePartitioning(nz, R, &nz_part);
    
    global_part = jxf_CTAlloc(JXF_Int, P*Q*R+1);
    global_part[0] = 0;
    cnt = 1;
    for (iz = 0; iz < R; iz ++)
    {
        nz_size = nz_part[iz+1] - nz_part[iz];
        for (iy = 0; iy < Q; iy ++)
        {
            ny_size = ny_part[iy+1] - ny_part[iy];
            for (ix = 0; ix < P; ix ++)
            {
                nx_size = nx_part[ix+1] - nx_part[ix];
                global_part[cnt] = global_part[cnt-1];
                global_part[cnt++] += nx_size * ny_size * nz_size;
            }
        }
    }
    
    nx_local = nx_part[p+1] - nx_part[p];
    ny_local = ny_part[q+1] - ny_part[q];
    nz_local = nz_part[r+1] - nz_part[r];
    
    my_id = r * (P * Q) + q * P + p;
    num_procs = P * Q * R;
    local_num_rows = nx_local * ny_local * nz_local;
    
    diag_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    offd_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    P_busy = jxf_min(nx, P);
    Q_busy = jxf_min(ny, Q);
    R_busy = jxf_min(nz, R);
    
    num_cols_offd = 0;
    if (p)
    {
        num_cols_offd += ny_local * nz_local;
    }
    if (p < P_busy-1)
    {
        num_cols_offd += ny_local * nz_local;
    }
    if (q)
    {
        num_cols_offd += nx_local * nz_local;
    }
    if (q < Q_busy-1)
    {
        num_cols_offd += nx_local * nz_local;
    }
    if (r)
    {
        num_cols_offd += nx_local * ny_local;
    }
    if (r < R_busy-1)
    {
        num_cols_offd += nx_local * ny_local;
    }
    if (!local_num_rows)
    {
        num_cols_offd = 0;
    }
    
    col_map_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
    cnt = 1;
    o_cnt = 1;
    diag_i[0] = 0;
    offd_i[0] = 0;
    for (iz = nz_part[r]; iz < nz_part[r+1]; iz ++)
    {
        for (iy = ny_part[q]; iy < ny_part[q+1]; iy ++)
        {
            for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
            {
                diag_i[cnt] = diag_i[cnt-1];
                offd_i[o_cnt] = offd_i[o_cnt-1];
                diag_i[cnt] ++;
                if (iz > nz_part[r])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iz)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iy > ny_part[q])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iy)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix > nx_part[p])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix+1 < nx_part[p+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iy+1 < ny_part[q+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iy+1 < ny)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iz+1 < nz_part[r+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iz+1 < nz)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                cnt ++;
                o_cnt ++;
            }
        }
    }
    diag_j = jxf_CTAlloc(JXF_Int, diag_i[local_num_rows]);
    diag_data = jxf_CTAlloc(JXF_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jxf_CTAlloc(JXF_Int, offd_i[local_num_rows]);
        offd_data = jxf_CTAlloc(JXF_Real, offd_i[local_num_rows]);
    }
    row_index = 0;
    cnt = 0;
    o_cnt = 0;
    for (iz = nz_part[r]; iz < nz_part[r+1]; iz ++)
    {
        for (iy = ny_part[q]; iy < ny_part[q+1]; iy ++)
        {
            for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
            {
                diag_j[cnt] = row_index;
                diag_data[cnt++] = value[0];
                if (iz > nz_part[r])
                {
                    diag_j[cnt] = row_index - nx_local * ny_local;
                    diag_data[cnt++] = value[3];
                }
                else
                {
                    if (iz)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy, iz-1, p, q, r-1, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[3];
                    }
                }
                if (iy > ny_part[q])
                {
                    diag_j[cnt] = row_index - nx_local;
                    diag_data[cnt++] = value[2];
                }
                else
                {
                    if (iy)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy-1, iz, p, q-1, r, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[2];
                    }
                }
                if (ix > nx_part[p])
                {
                    diag_j[cnt] = row_index - 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix)
                    {
                        offd_j[o_cnt] = jxf_map(ix-1, iy, iz, p-1, q, r, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                if (ix+1 < nx_part[p+1])
                {
                    diag_j[cnt] = row_index + 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jxf_map(ix+1, iy, iz, p+1, q, r, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                if (iy+1 < ny_part[q+1])
                {
                    diag_j[cnt] = row_index + nx_local;
                    diag_data[cnt++] = value[2];
                }
                else
                {
                    if (iy+1 < ny)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy+1, iz, p, q+1, r, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[2];
                    }
                }
                if (iz+1 < nz_part[r+1])
                {
                    diag_j[cnt] = row_index + nx_local * ny_local;
                    diag_data[cnt++] = value[3];
                }
                else
                {
                    if (iz+1 < nz)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy, iz+1, p, q, r+1, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[3];
                    }
                }
                row_index ++;
            }
        }
    }
    if (num_procs > 1)
    {
        for (i = 0; i < num_cols_offd; i ++)
        {
            col_map_offd[i] = offd_j[i];
        }
        jxf_qsort0(col_map_offd, 0, num_cols_offd-1);
        for (i = 0; i < num_cols_offd; i ++)
        {
            for (j = 0; j < num_cols_offd; j ++)
            {
                if (offd_j[i] == col_map_offd[j])
                {
                    offd_j[i] = j;
                    break;
                }
            }
        }
    }
#ifdef JXF_NO_GLOBAL_PARTITION
/* ideally we would use less storage earlier in this function, but this is fine for testing */
    {
        JXF_Int tmp1, tmp2;
        tmp1 = global_part[my_id];
        tmp2 = global_part[my_id+1];
        jxf_TFree(global_part);
        global_part = jxf_CTAlloc(JXF_Int, 2);
        global_part[0] = tmp1;
        global_part[1] = tmp2;
    }
#endif
    A = jxf_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                      global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jxf_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jxf_ParCSRMatrixDiag(A);
    jxf_CSRMatrixI(diag) = diag_i;
    jxf_CSRMatrixJ(diag) = diag_j;
    jxf_CSRMatrixData(diag) = diag_data;
    offd = jxf_ParCSRMatrixOffd(A);
    jxf_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jxf_CSRMatrixJ(offd) = offd_j;
        jxf_CSRMatrixData(offd) = offd_data;
    }
    jxf_TFree(nx_part);
    jxf_TFree(ny_part);
    jxf_TFree(nz_part);
    
    return A;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_GenerateParLaplacian2d9pt
 */
jxf_ParCSRMatrix *
jxf_GenerateParLaplacian2d9pt( MPI_Comm comm,
                              JXF_Int nx,
                              JXF_Int ny,
                              JXF_Int P,
                              JXF_Int Q,
                              JXF_Int p,
                              JXF_Int q,
                              JXF_Real *value )
{
    jxf_ParCSRMatrix *A = NULL;
    
    jxf_CSRMatrix *diag = NULL;
    jxf_CSRMatrix *offd = NULL;
    
    JXF_Int *diag_i = NULL;
    JXF_Int *diag_j = NULL;
    JXF_Real *diag_data = NULL;
    
    JXF_Int *offd_i = NULL;
    JXF_Int *offd_j = NULL;
    JXF_Real *offd_data = NULL;
    
    JXF_Int *work = NULL;
    JXF_Int *nx_part = NULL;
    JXF_Int *ny_part = NULL;
    JXF_Int *global_part = NULL;
    JXF_Int *col_map_offd = NULL;
    
    JXF_Int i, j, ix, iy, cnt, o_cnt;
    JXF_Int num_cols_offd, grid_size;
    JXF_Int local_num_rows, row_index;
    JXF_Int num_procs, my_id, P_busy, Q_busy;
    JXF_Int nx_local, ny_local, nx_size, ny_size;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny;
    
    jxf_GeneratePartitioning(nx, P, &nx_part);
    jxf_GeneratePartitioning(ny, Q, &ny_part);
    
    global_part = jxf_CTAlloc(JXF_Int, P*Q+1);
    global_part[0] = 0;
    cnt = 1;
    for (iy = 0; iy < Q; iy ++)
    {
        ny_size = ny_part[iy+1] - ny_part[iy];
        for (ix = 0; ix < P; ix ++)
        {
            nx_size = nx_part[ix+1] - nx_part[ix];
            global_part[cnt] = global_part[cnt-1];
            global_part[cnt++] += nx_size * ny_size;
        }
    }
    nx_local = nx_part[p+1] - nx_part[p];
    ny_local = ny_part[q+1] - ny_part[q];
    
    my_id = q * P + p;
    num_procs = P * Q;
    local_num_rows = nx_local * ny_local;
    
    diag_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    offd_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    P_busy = jxf_min(nx, P);
    Q_busy = jxf_min(ny, Q);
    
    num_cols_offd = 0;
    if (p)
    {
        num_cols_offd += ny_local;
    }
    if (p < P_busy-1)
    {
        num_cols_offd += ny_local;
    }
    if (q)
    {
        num_cols_offd += nx_local;
    }
    if (q < Q_busy-1)
    {
        num_cols_offd += nx_local;
    }
    if (p && q)
    {
        num_cols_offd ++;
    }
    if (p && q < Q_busy-1)
    {
        num_cols_offd ++;
    }
    if (p < P_busy-1 && q)
    {
        num_cols_offd ++;
    }
    if (p < P_busy-1 && q < Q_busy-1)
    {
        num_cols_offd ++;
    }
    if (!local_num_rows)
    {
        num_cols_offd = 0;
    }
    col_map_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
    cnt = 0;
    o_cnt = 0;
    diag_i[0] = 0;
    offd_i[0] = 0;
    for (iy = ny_part[q];  iy < ny_part[q+1]; iy ++)
    {
        for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
        {
            cnt ++;
            o_cnt ++;
            diag_i[cnt] = diag_i[cnt-1];
            offd_i[o_cnt] = offd_i[o_cnt-1];
            diag_i[cnt] ++;
            if (iy > ny_part[q])
            {
                diag_i[cnt] ++;
                if (ix > nx_part[p])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix < nx_part[p+1]-1)
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
            }
            else
            {
                if (iy)
                {
                    offd_i[o_cnt] ++;
                    if (ix > nx_part[p])
                    {
                        offd_i[o_cnt] ++;
                    }
                    else if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                    if (ix < nx_part[p+1]-1)
                    {
                        offd_i[o_cnt] ++;
                    }
                    else if (ix < nx-1)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
            }
            if (ix > nx_part[p])
            {
                diag_i[cnt] ++;
            }
            else
            {
                if (ix)
                {
                    offd_i[o_cnt] ++;
                }
            }
            if (ix+1 < nx_part[p+1])
            {
                diag_i[cnt] ++;
            }
            else
            {
                if (ix+1 < nx)
                {
                    offd_i[o_cnt] ++;
                }
            }
            if (iy+1 < ny_part[q+1])
            {
                diag_i[cnt] ++;
                if (ix > nx_part[p])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix < nx_part[p+1]-1)
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
            }
            else
            {
                if (iy+1 < ny)
                {
                    offd_i[o_cnt] ++;
                    if (ix > nx_part[p])
                    {
                        offd_i[o_cnt] ++;
                    }
                    else if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                    if (ix < nx_part[p+1]-1)
                    {
                        offd_i[o_cnt] ++;
                    }
                    else if (ix < nx-1)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
            }
        }
    }
    diag_j = jxf_CTAlloc(JXF_Int, diag_i[local_num_rows]);
    diag_data = jxf_CTAlloc(JXF_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jxf_CTAlloc(JXF_Int, offd_i[local_num_rows]);
        offd_data = jxf_CTAlloc(JXF_Real, offd_i[local_num_rows]);
    }
    row_index = 0;
    cnt = 0;
    o_cnt = 0;
    for (iy = ny_part[q];  iy < ny_part[q+1]; iy ++)
    {
        for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
        {
            diag_j[cnt] = row_index;
            diag_data[cnt++] = value[0];
            if (iy > ny_part[q])
            {
                if (ix > nx_part[p])
                {
                    diag_j[cnt] = row_index - nx_local - 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix)
                    {
                        offd_j[o_cnt] = jxf_map2(ix-1, iy-1, p-1, q, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                diag_j[cnt] = row_index - nx_local;
                diag_data[cnt++] = value[1];
                if (ix < nx_part[p+1]-1)
                {
                    diag_j[cnt] = row_index - nx_local + 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jxf_map2(ix+1, iy-1, p+1, q, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
            }
            else
            {
                if (iy)
                {
                    if (ix > nx_part[p])
                    {
                        offd_j[o_cnt] = jxf_map2(ix-1, iy-1, p, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix)
                    {
                        offd_j[o_cnt] = jxf_map2(ix-1, iy-1, p-1, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    offd_j[o_cnt] = jxf_map2(ix, iy-1, p, q-1, P, Q, nx_part, ny_part, global_part);
                    offd_data[o_cnt++] = value[1];
                    if (ix < nx_part[p+1]-1)
                    {
                        offd_j[o_cnt] = jxf_map2(ix+1, iy-1, p, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jxf_map2(ix+1, iy-1, p+1, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
            }
            if (ix > nx_part[p])
            {
                diag_j[cnt] = row_index - 1;
                diag_data[cnt++] = value[1];
            }
            else
            {
                if (ix)
                {
                    offd_j[o_cnt] = jxf_map2(ix-1, iy, p-1, q, P, Q, nx_part, ny_part, global_part);
                    offd_data[o_cnt++] = value[1];
                }
            }
            if (ix+1 < nx_part[p+1])
            {
                diag_j[cnt] = row_index + 1;
                diag_data[cnt++] = value[1];
            }
            else
            {
                if (ix+1 < nx)
                {
                    offd_j[o_cnt] = jxf_map2(ix+1, iy, p+1, q, P, Q, nx_part, ny_part, global_part);
                    offd_data[o_cnt++] = value[1];
                }
            }
            if (iy+1 < ny_part[q+1])
            {
                if (ix > nx_part[p])
                {
                    diag_j[cnt] = row_index + nx_local - 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix)
                    {
                        offd_j[o_cnt] = jxf_map2(ix-1, iy+1, p-1, q, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                diag_j[cnt] = row_index + nx_local;
                diag_data[cnt++] = value[1];
                if (ix < nx_part[p+1]-1)
                {
                    diag_j[cnt] = row_index + nx_local + 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jxf_map2(ix+1, iy+1, p+1, q, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
            }
            else
            {
                if (iy+1 < ny)
                {
                    if (ix > nx_part[p])
                    {
                        offd_j[o_cnt] = jxf_map2(ix-1, iy+1, p, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix)
                    {
                        offd_j[o_cnt] = jxf_map2(ix-1, iy+1, p-1, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    offd_j[o_cnt] = jxf_map2(ix, iy+1, p, q+1, P, Q, nx_part, ny_part, global_part);
                    offd_data[o_cnt++] = value[1];
                    if (ix < nx_part[p+1]-1)
                    {
                        offd_j[o_cnt] = jxf_map2(ix+1, iy+1, p, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix < nx-1)
                    {
                        offd_j[o_cnt] = jxf_map2(ix+1, iy+1, p+1, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
            }
            row_index ++;
        }
    }
    if (num_procs > 1)
    {
        work = jxf_CTAlloc(JXF_Int, o_cnt);
        for (i = 0; i < o_cnt; i ++)
        {
            work[i] = offd_j[i];
        }
        jxf_qsort0(work, 0, o_cnt-1);
        col_map_offd[0] = work[0];
        cnt = 0;
        for (i = 0; i < o_cnt; i ++)
        {
            if (work[i] > col_map_offd[cnt])
            {
                cnt ++;
                col_map_offd[cnt] = work[i];
            }
        }
        for (i = 0; i < o_cnt; i ++)
        {
            for (j = 0; j < num_cols_offd; j ++)
            {
                if (offd_j[i] == col_map_offd[j])
                {
                    offd_j[i] = j;
                    break;
                }
            }
        }
        jxf_TFree(work);
    }
#ifdef JXF_NO_GLOBAL_PARTITION
/* ideally we would use less storage earlier in this function, but this is fine for testing */
    {
        JXF_Int tmp1, tmp2;
        tmp1 = global_part[my_id];
        tmp2 = global_part[my_id+1];
        jxf_TFree(global_part);
        global_part = jxf_CTAlloc(JXF_Int, 2);
        global_part[0] = tmp1;
        global_part[1] = tmp2;
    }
#endif
    A = jxf_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                 global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jxf_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jxf_ParCSRMatrixDiag(A);
    jxf_CSRMatrixI(diag) = diag_i;
    jxf_CSRMatrixJ(diag) = diag_j;
    jxf_CSRMatrixData(diag) = diag_data;
    offd = jxf_ParCSRMatrixOffd(A);
    jxf_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jxf_CSRMatrixJ(offd) = offd_j;
        jxf_CSRMatrixData(offd) = offd_data;
    }
    jxf_TFree(nx_part);
    jxf_TFree(ny_part);
    
    return A;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_GenerateParLaplacian2d9pt
 */
jxf_ParCSRMatrix *
jxf_GenerateParLaplacian3d27pt( MPI_Comm comm,
                               JXF_Int nx,
                               JXF_Int ny,
                               JXF_Int nz,
                               JXF_Int P,
                               JXF_Int Q,
                               JXF_Int R,
                               JXF_Int p,
                               JXF_Int q,
                               JXF_Int r,
                               JXF_Real *value )
{
    jxf_ParCSRMatrix *A = NULL;
    
    jxf_CSRMatrix *diag = NULL;
    jxf_CSRMatrix *offd = NULL;
    
    JXF_Int *diag_i = NULL;
    JXF_Int *diag_j = NULL;
    JXF_Real *diag_data = NULL;
    
    JXF_Int *offd_i = NULL;
    JXF_Int *offd_j = NULL;
    JXF_Real *offd_data = NULL;
    
    JXF_Int *work = NULL;
    JXF_Int *nx_part = NULL;
    JXF_Int *ny_part = NULL;
    JXF_Int *nz_part = NULL;
    JXF_Int *global_part = NULL;
    JXF_Int *col_map_offd = NULL;
    
    JXF_Int i, j, ix, iy, iz;
    JXF_Int num_cols_offd, nxy, grid_size;
    JXF_Int cnt, o_cnt, local_num_rows, row_index;
    JXF_Int num_procs, my_id, P_busy, Q_busy, R_busy;
    JXF_Int nx_local, ny_local, nz_local, nx_size, ny_size, nz_size;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny * nz;
    
    jxf_GeneratePartitioning(nx, P, &nx_part);
    jxf_GeneratePartitioning(ny, Q, &ny_part);
    jxf_GeneratePartitioning(nz, R, &nz_part);
    
    global_part = jxf_CTAlloc(JXF_Int, P*Q*R+1);
    global_part[0] = 0;
    cnt = 1;
    for (iz = 0; iz < R; iz ++)
    {
        nz_size = nz_part[iz+1] - nz_part[iz];
        for (iy = 0; iy < Q; iy ++)
        {
            ny_size = (ny_part[iy+1] - ny_part[iy]) * nz_size;
            for (ix = 0; ix < P; ix ++)
            {
                nx_size = nx_part[ix+1] - nx_part[ix];
                global_part[cnt] = global_part[cnt-1];
                global_part[cnt++] += nx_size * ny_size;
            }
        }
    }
    nx_local = nx_part[p+1] - nx_part[p];
    ny_local = ny_part[q+1] - ny_part[q];
    nz_local = nz_part[r+1] - nz_part[r];
    
    my_id = (r * Q + q) * P + p;
    num_procs = P * Q * R;
    
    local_num_rows = nx_local * ny_local * nz_local;
    diag_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    offd_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    
    P_busy = jxf_min(nx, P);
    Q_busy = jxf_min(ny, Q);
    R_busy = jxf_min(nz, R);
    
    num_cols_offd = 0;
    if (p)
    {
        num_cols_offd += ny_local * nz_local;
    }
    if (p < P_busy-1)
    {
        num_cols_offd += ny_local * nz_local;
    }
    if (q)
    {
        num_cols_offd += nx_local * nz_local;
    }
    if (q < Q_busy-1)
    {
        num_cols_offd += nx_local * nz_local;
    }
    if (r)
    {
        num_cols_offd += nx_local * ny_local;
    }
    if (r < R_busy-1)
    {
        num_cols_offd += nx_local * ny_local;
    }
    if (p && q)
    {
        num_cols_offd += nz_local;
    }
    if (p && q < Q_busy-1)
    {
        num_cols_offd += nz_local;
    }
    if (p < P_busy-1 && q)
    {
        num_cols_offd += nz_local;
    }
    if (p < P_busy-1 && q < Q_busy-1)
    {
        num_cols_offd += nz_local;
    }
    if (p && r)
    {
        num_cols_offd += ny_local;
    }
    if (p && r < R_busy-1)
    {
        num_cols_offd += ny_local;
    }
    if (p < P_busy-1 && r)
    {
        num_cols_offd += ny_local;
    }
    if (p < P_busy-1 && r < R_busy-1)
    {
        num_cols_offd += ny_local;
    }
    if (q && r)
    {
        num_cols_offd += nx_local;
    }
    if (q && r < R_busy-1)
    {
        num_cols_offd += nx_local;
    }
    if (q < Q_busy-1 && r)
    {
        num_cols_offd += nx_local;
    }
    if (q < Q_busy-1 && r < R_busy-1)
    {
        num_cols_offd += nx_local;
    }
    if (p && q && r)
    {
        num_cols_offd ++;
    }
    if (p && q && r < R_busy-1)
    {
        num_cols_offd ++;
    }
    if (p && q < Q_busy-1 && r)
    {
        num_cols_offd ++;
    }
    if (p && q < Q_busy-1 && r < R_busy-1)
    {
        num_cols_offd ++;
    }
    if (p < P_busy-1 && q && r)
    {
        num_cols_offd ++;
    }
    if (p < P_busy-1 && q && r < R_busy-1)
    {
        num_cols_offd ++;
    }
    if (p < P_busy-1 && q < Q_busy-1 && r)
    {
        num_cols_offd ++;
    }
    if (p < P_busy-1 && q < Q_busy-1 && r < R_busy-1)
    {
        num_cols_offd ++;
    }
    if (!local_num_rows)
    {
        num_cols_offd = 0;
    }
    col_map_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
    cnt = 0;
    o_cnt = 0;
    diag_i[0] = 0;
    offd_i[0] = 0;
    for (iz = nz_part[r];  iz < nz_part[r+1]; iz ++)
    {
        for (iy = ny_part[q];  iy < ny_part[q+1]; iy ++)
        {
            for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
            {
                cnt ++;
                o_cnt ++;
                diag_i[cnt] = diag_i[cnt-1];
                offd_i[o_cnt] = offd_i[o_cnt-1];
                diag_i[cnt] ++;
                if (iz > nz_part[r])
                {
                    diag_i[cnt] ++;
                    if (iy > ny_part[q])
                    {
                        diag_i[cnt] ++;
                        if (ix > nx_part[p])
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                    else
                    {
                        if (iy)
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix < nx-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                    if (ix > nx_part[p])
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                    if (ix+1 < nx_part[p+1])
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                    if (iy+1 < ny_part[q+1])
                    {
                        diag_i[cnt] ++;
                        if (ix > nx_part[p])
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                    else
                    {
                        if (iy+1 < ny)
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix < nx-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                }
                else
                {
                    if (iz)
                    {
                        offd_i[o_cnt] ++;
                        if (iy > ny_part[q])
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                        else
                        {
                            if (iy)
                            {
                                offd_i[o_cnt] ++;
                                if (ix > nx_part[p])
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix < nx-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                        if (ix > nx_part[p])
                        {
                            offd_i[o_cnt] ++;
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (ix+1 < nx_part[p+1])
                        {
                            offd_i[o_cnt] ++;
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (iy+1 < ny_part[q+1])
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                        else
                        {
                            if (iy+1 < ny)
                            {
                                offd_i[o_cnt] ++;
                                if (ix > nx_part[p])
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix < nx-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                    }
                }
                if (iy > ny_part[q])
                {
                    diag_i[cnt] ++;
                    if (ix > nx_part[p])
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                    if (ix < nx_part[p+1]-1)
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                }
                else
                {
                    if (iy)
                    {
                        offd_i[o_cnt] ++;
                        if (ix > nx_part[p])
                        {
                            offd_i[o_cnt] ++;
                        }
                        else if (ix)
                        {
                            offd_i[o_cnt] ++;
                        }
                        if (ix < nx_part[p+1]-1)
                        {
                            offd_i[o_cnt] ++;
                        }
                        else if (ix < nx-1)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                }
                if (ix > nx_part[p])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix+1 < nx_part[p+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iy+1 < ny_part[q+1])
                {
                    diag_i[cnt] ++;
                    if (ix > nx_part[p])
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                    if (ix < nx_part[p+1]-1)
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                }
                else
                {
                    if (iy+1 < ny)
                    {
                        offd_i[o_cnt] ++;
                        if (ix > nx_part[p])
                        {
                            offd_i[o_cnt] ++;
                        }
                        else if (ix)
                        {
                            offd_i[o_cnt] ++;
                        }
                        if (ix < nx_part[p+1]-1)
                        {
                            offd_i[o_cnt] ++;
                        }
                        else if (ix < nx-1)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                }
                if (iz+1 < nz_part[r+1])
                {
                    diag_i[cnt] ++;
                    if (iy > ny_part[q])
                    {
                        diag_i[cnt] ++;
                        if (ix > nx_part[p])
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                    else
                    {
                        if (iy)
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix < nx-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                    if (ix > nx_part[p])
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                    if (ix+1 < nx_part[p+1])
                    {
                        diag_i[cnt] ++;
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_i[o_cnt] ++;
                        }
                    }
                    if (iy+1 < ny_part[q+1])
                    {
                        diag_i[cnt] ++;
                        if (ix > nx_part[p])
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_i[cnt] ++;
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                    else
                    {
                        if (iy+1 < ny)
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else if (ix < nx-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                    }
                }
                else
                {
                    if (iz+1 < nz)
                    {
                        offd_i[o_cnt] ++;
                        if (iy > ny_part[q])
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                        else
                        {
                            if (iy)
                            {
                                offd_i[o_cnt] ++;
                                if (ix > nx_part[p])
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix < nx-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                        if (ix > nx_part[p])
                        {
                            offd_i[o_cnt] ++;
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (ix+1 < nx_part[p+1])
                        {
                            offd_i[o_cnt] ++;
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_i[o_cnt] ++;
                            }
                        }
                        if (iy+1 < ny_part[q+1])
                        {
                            offd_i[o_cnt] ++;
                            if (ix > nx_part[p])
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_i[o_cnt] ++;
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                        else
                        {
                            if (iy+1 < ny)
                            {
                                offd_i[o_cnt] ++;
                                if (ix > nx_part[p])
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                                else if (ix < nx-1)
                                {
                                    offd_i[o_cnt] ++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    diag_j = jxf_CTAlloc(JXF_Int, diag_i[local_num_rows]);
    diag_data = jxf_CTAlloc(JXF_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jxf_CTAlloc(JXF_Int, offd_i[local_num_rows]);
        offd_data = jxf_CTAlloc(JXF_Real, offd_i[local_num_rows]);
    }
    nxy = nx_local * ny_local;
    row_index = 0;
    cnt = 0;
    o_cnt = 0;
    for (iz = nz_part[r];  iz < nz_part[r+1]; iz ++)
    {
        for (iy = ny_part[q];  iy < ny_part[q+1]; iy ++)
        {
            for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
            {
                diag_j[cnt] = row_index;
                diag_data[cnt++] = value[0];
                if (iz > nz_part[r])
                {
                    if (iy > ny_part[q])
                    {
                        if (ix > nx_part[p])
                        {
                            diag_j[cnt] = row_index - nxy - nx_local - 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p-1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        diag_j[cnt] = row_index - nxy - nx_local;
                        diag_data[cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_j[cnt] = row_index - nxy - nx_local + 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p+1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                    else
                    {
                        if (iy)
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p-1, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy-1, iz-1, p, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p+1, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                    if (ix > nx_part[p])
                    {
                        diag_j[cnt] = row_index - nxy - 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy, iz-1, p-1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                    diag_j[cnt] = row_index - nxy;
                    diag_data[cnt++] = value[1];
                    if (ix+1 < nx_part[p+1])
                    {
                        diag_j[cnt] = row_index - nxy + 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy, iz-1, p+1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                    if (iy+1 < ny_part[q+1])
                    {
                        if (ix > nx_part[p])
                        {
                            diag_j[cnt] = row_index - nxy + nx_local - 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p-1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        diag_j[cnt] = row_index - nxy + nx_local;
                        diag_data[cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_j[cnt] = row_index - nxy + nx_local + 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p+1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                    else
                    {
                        if (iy+1 < ny)
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p-1, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy+1, iz-1, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p+1, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                }
                else
                {
                    if (iz)
                    {
                        if (iy > ny_part[q])
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p-1, q, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy-1, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p+1, q, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        else
                        {
                            if (iy)
                            {
                                if (ix > nx_part[p])
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz-1, p-1, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jxf_map3(ix, iy-1, iz-1, p, q-1, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz-1, p+1, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        if (ix > nx_part[p])
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy, iz-1, p-1, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        offd_j[o_cnt] = jxf_map3(ix, iy, iz-1, p, q, r-1,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix+1 < nx_part[p+1])
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy, iz-1, p+1, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        if (iy+1 < ny_part[q+1])
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p-1, q, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy+1, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p+1, q, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        else
                        {
                            if (iy+1 < ny)
                            {
                                if (ix > nx_part[p])
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz-1, p-1, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jxf_map3(ix, iy+1, iz-1, p, q+1, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz-1, p+1, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                    }
                }
                if (iy > ny_part[q])
                {
                    if (ix > nx_part[p])
                    {
                        diag_j[cnt] = row_index - nx_local - 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz, p-1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                    diag_j[cnt] = row_index - nx_local;
                    diag_data[cnt++] = value[1];
                    if (ix < nx_part[p+1]-1)
                    {
                        diag_j[cnt] = row_index - nx_local + 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz, p+1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                }
                else
                {
                    if (iy)
                    {
                        if (ix > nx_part[p])
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz, p, q-1, r,
                                               P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix)
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz, p-1, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        offd_j[o_cnt] = jxf_map3(ix, iy-1, iz, p, q-1, r,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz, p, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix < nx-1)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz, p+1, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                }
                if (ix > nx_part[p])
                {
                    diag_j[cnt] = row_index - 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix)
                    {
                        offd_j[o_cnt] = jxf_map3(ix-1, iy, iz, p-1, q, r,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                if (ix+1 < nx_part[p+1])
                {
                    diag_j[cnt] = row_index + 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jxf_map3(ix+1, iy, iz, p+1, q, r,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                if (iy+1 < ny_part[q+1])
                {
                    if (ix > nx_part[p])
                    {
                        diag_j[cnt] = row_index + nx_local - 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz, p-1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                    diag_j[cnt] = row_index + nx_local;
                    diag_data[cnt++] = value[1];
                    if (ix < nx_part[p+1]-1)
                    {
                        diag_j[cnt] = row_index + nx_local + 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz, p+1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                }
                else
                {
                    if (iy+1 < ny)
                    {
                        if (ix > nx_part[p])
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix)
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz, p-1, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        offd_j[o_cnt] = jxf_map3(ix, iy+1, iz, p, q+1, r,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix < nx-1)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz, p+1, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                }
                if (iz+1 < nz_part[r+1])
                {
                    if (iy > ny_part[q])
                    {
                        if (ix > nx_part[p])
                        {
                            diag_j[cnt] = row_index + nxy - nx_local - 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p-1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        diag_j[cnt] = row_index + nxy - nx_local;
                        diag_data[cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_j[cnt] = row_index + nxy - nx_local + 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p+1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                    else
                    {
                        if (iy)
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p-1, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy-1, iz+1, p, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p+1, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                    if (ix > nx_part[p])
                    {
                        diag_j[cnt] = row_index + nxy - 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix)
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy, iz+1, p-1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                    diag_j[cnt] = row_index + nxy;
                    diag_data[cnt++] = value[1];
                    if (ix+1 < nx_part[p+1])
                    {
                        diag_j[cnt] = row_index + nxy + 1;
                        diag_data[cnt++] = value[1];
                    }
                    else
                    {
                        if (ix+1 < nx)
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy, iz+1, p+1, q, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                    }
                    if (iy+1 < ny_part[q+1])
                    {
                        if (ix > nx_part[p])
                        {
                            diag_j[cnt] = row_index + nxy + nx_local - 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p-1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        diag_j[cnt] = row_index + nxy + nx_local;
                        diag_data[cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            diag_j[cnt] = row_index + nxy + nx_local + 1;
                            diag_data[cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p+1, q, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                    else
                    {
                        if (iy+1 < ny)
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p-1, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy+1, iz+1, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p+1, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                    }
                }
                else
                {
                    if (iz+1 < nz)
                    {
                        if (iy > ny_part[q])
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p-1, q, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy-1, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p+1, q, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        else
                        {
                            if (iy)
                            {
                                if (ix > nx_part[p])
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy-1, iz+1, p-1, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jxf_map3(ix, iy-1, iz+1, p, q-1, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy-1, iz+1, p+1, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        if (ix > nx_part[p])
                        {
                            offd_j[o_cnt] = jxf_map3(ix-1, iy, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy, iz+1, p-1, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        offd_j[o_cnt] = jxf_map3(ix, iy, iz+1, p, q, r+1,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix+1 < nx_part[p+1])
                        {
                            offd_j[o_cnt] = jxf_map3(ix+1, iy, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy, iz+1, p+1, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        if (iy+1 < ny_part[q+1])
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p-1, q, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jxf_map3(ix, iy+1, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p+1, q, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        else
                        {
                            if (iy+1 < ny)
                            {
                                if (ix > nx_part[p])
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix-1, iy+1, iz+1, p-1, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jxf_map3(ix, iy+1, iz+1, p, q+1, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jxf_map3(ix+1, iy+1, iz+1, p+1, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                    }
                }
                row_index ++;
            }
        }
    }
    if (num_procs > 1)
    {
        work = jxf_CTAlloc(JXF_Int, o_cnt);
        for (i = 0; i < o_cnt; i ++)
        {
            work[i] = offd_j[i];
        }
        jxf_qsort0(work, 0, o_cnt-1);
        col_map_offd[0] = work[0];
        cnt = 0;
        for (i = 0; i < o_cnt; i ++)
        {
            if (work[i] > col_map_offd[cnt])
            {
                cnt ++;
                col_map_offd[cnt] = work[i];
            }
        }
        for (i = 0; i < o_cnt; i ++)
        {
            for (j = 0; j < num_cols_offd; j ++)
            {
                if (offd_j[i] == col_map_offd[j])
                {
                    offd_j[i] = j;
                    break;
                }
            }
        }
        jxf_TFree(work);
    }
#ifdef JXF_NO_GLOBAL_PARTITION
/* ideally we would use less storage earlier in this function, but this is fine for testing */
    {
        JXF_Int tmp1, tmp2;
        tmp1 = global_part[my_id];
        tmp2 = global_part[my_id+1];
        jxf_TFree(global_part);
        global_part = jxf_CTAlloc(JXF_Int, 2);
        global_part[0] = tmp1;
        global_part[1] = tmp2;
    }
#endif
    A = jxf_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                   global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jxf_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jxf_ParCSRMatrixDiag(A);
    jxf_CSRMatrixI(diag) = diag_i;
    jxf_CSRMatrixJ(diag) = diag_j;
    jxf_CSRMatrixData(diag) = diag_data;
    offd = jxf_ParCSRMatrixOffd(A);
    jxf_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jxf_CSRMatrixJ(offd) = offd_j;
        jxf_CSRMatrixData(offd) = offd_data;
    }
    jxf_TFree(nx_part);
    jxf_TFree(ny_part);
    jxf_TFree(nz_part);
    
    return A;
}

/*!
 * \fn jxf_ParCSRMatrix *jxf_GenerateParConvecDiff
 */
jxf_ParCSRMatrix *
jxf_GenerateParConvecDiff( MPI_Comm comm,
                          JXF_Int nx,
                          JXF_Int ny,
                          JXF_Int nz,
                          JXF_Int P,
                          JXF_Int Q,
                          JXF_Int R,
                          JXF_Int p,
                          JXF_Int q,
                          JXF_Int r,
                          JXF_Real  *value )
{
    jxf_ParCSRMatrix *A = NULL;
    
    jxf_CSRMatrix *diag = NULL;
    jxf_CSRMatrix *offd = NULL;
    
    JXF_Int *diag_i = NULL;
    JXF_Int *diag_j = NULL;
    JXF_Real *diag_data = NULL;
    
    JXF_Int *offd_i = NULL;
    JXF_Int *offd_j = NULL;
    JXF_Real *offd_data = NULL;
    
    JXF_Int *nx_part = NULL;
    JXF_Int *ny_part = NULL;
    JXF_Int *nz_part = NULL;
    JXF_Int *global_part = NULL;
    JXF_Int *col_map_offd = NULL;
    
    JXF_Int i, j, ix, iy, iz;
    JXF_Int num_cols_offd, grid_size;
    JXF_Int row_index, cnt, o_cnt, local_num_rows;
    JXF_Int num_procs, my_id, P_busy, Q_busy, R_busy;
    JXF_Int nx_local, ny_local, nz_local, nx_size, ny_size, nz_size;
    
    jxf_MPI_Comm_size(comm, &num_procs);
    jxf_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny * nz;
    
    jxf_GeneratePartitioning(nx, P, &nx_part);
    jxf_GeneratePartitioning(ny, Q, &ny_part);
    jxf_GeneratePartitioning(nz, R, &nz_part);
    
    global_part = jxf_CTAlloc(JXF_Int, P*Q*R+1);
    
    global_part[0] = 0;
    cnt = 1;
    for (iz = 0; iz < R; iz ++)
    {
        nz_size = nz_part[iz+1] - nz_part[iz];
        for (iy = 0; iy < Q; iy ++)
        {
            ny_size = ny_part[iy+1] - ny_part[iy];
            for (ix = 0; ix < P; ix ++)
            {
                nx_size = nx_part[ix+1] - nx_part[ix];
                global_part[cnt] = global_part[cnt-1];
                global_part[cnt++] += nx_size * ny_size * nz_size;
            }
        }
    }
    nx_local = nx_part[p+1] - nx_part[p];
    ny_local = ny_part[q+1] - ny_part[q];
    nz_local = nz_part[r+1] - nz_part[r];
    
    my_id = r * (P * Q) + q * P + p;
    num_procs = P * Q * R;
    
    local_num_rows = nx_local * ny_local * nz_local;
    diag_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    offd_i = jxf_CTAlloc(JXF_Int, local_num_rows+1);
    
    P_busy = jxf_min(nx, P);
    Q_busy = jxf_min(ny, Q);
    R_busy = jxf_min(nz, R);
    
    num_cols_offd = 0;
    if (p)
    {
        num_cols_offd += ny_local * nz_local;
    }
    if (p < P_busy-1)
    {
        num_cols_offd += ny_local * nz_local;
    }
    if (q)
    {
        num_cols_offd += nx_local * nz_local;
    }
    if (q < Q_busy-1)
    {
        num_cols_offd += nx_local * nz_local;
    }
    if (r)
    {
        num_cols_offd += nx_local * ny_local;
    }
    if (r < R_busy-1)
    {
        num_cols_offd += nx_local * ny_local;
    }
    if (!local_num_rows)
    {
        num_cols_offd = 0;
    }
    col_map_offd = jxf_CTAlloc(JXF_Int, num_cols_offd);
    cnt = 1;
    o_cnt = 1;
    diag_i[0] = 0;
    offd_i[0] = 0;
    for (iz = nz_part[r]; iz < nz_part[r+1]; iz ++)
    {
        for (iy = ny_part[q];  iy < ny_part[q+1]; iy ++)
        {
            for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
            {
                diag_i[cnt] = diag_i[cnt-1];
                offd_i[o_cnt] = offd_i[o_cnt-1];
                diag_i[cnt] ++;
                if (iz > nz_part[r])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iz)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iy > ny_part[q])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iy)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix > nx_part[p])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (ix+1 < nx_part[p+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iy+1 < ny_part[q+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iy+1 < ny)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                if (iz+1 < nz_part[r+1])
                {
                    diag_i[cnt] ++;
                }
                else
                {
                    if (iz+1 < nz)
                    {
                        offd_i[o_cnt] ++;
                    }
                }
                cnt ++;
                o_cnt ++;
            }
        }
    }
    diag_j = jxf_CTAlloc(JXF_Int, diag_i[local_num_rows]);
    diag_data = jxf_CTAlloc(JXF_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jxf_CTAlloc(JXF_Int, offd_i[local_num_rows]);
        offd_data = jxf_CTAlloc(JXF_Real, offd_i[local_num_rows]);
    }
    row_index = 0;
    cnt = 0;
    o_cnt = 0;
    for (iz = nz_part[r]; iz < nz_part[r+1]; iz ++)
    {
        for (iy = ny_part[q];  iy < ny_part[q+1]; iy ++)
        {
            for (ix = nx_part[p]; ix < nx_part[p+1]; ix ++)
            {
                diag_j[cnt] = row_index;
                diag_data[cnt++] = value[0];
                if (iz > nz_part[r])
                {
                    diag_j[cnt] = row_index - nx_local * ny_local;
                    diag_data[cnt++] = value[3];
                }
                else
                {
                    if (iz)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy, iz-1, p, q, r-1, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[3];
                    }
                }
                if (iy > ny_part[q])
                {
                    diag_j[cnt] = row_index - nx_local;
                    diag_data[cnt++] = value[2];
                }
                else
                {
                    if (iy)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy-1, iz, p, q-1, r, P, Q, R,
                                                          nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[2];
                    }
                }
                if (ix > nx_part[p])
                {
                    diag_j[cnt] = row_index - 1;
                    diag_data[cnt++] = value[1];
                }
                else
                {
                    if (ix)
                    {
                        offd_j[o_cnt] = jxf_map(ix-1, iy, iz, p-1, q, r, P, Q, R,
                                                           nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
                if (ix+1 < nx_part[p+1])
                {
                    diag_j[cnt] = row_index + 1;
                    diag_data[cnt++] = value[4];
                }
                else
                {
                    if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jxf_map(ix+1, iy, iz, p+1, q, r, P, Q, R,
                                                          nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[4];
                    }
                }
                if (iy+1 < ny_part[q+1])
                {
                    diag_j[cnt] = row_index + nx_local;
                    diag_data[cnt++] = value[5];
                }
                else
                {
                    if (iy+1 < ny)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy+1, iz, p, q+1, r, P, Q, R,
                                                          nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[5];
                    }
                }
                if (iz+1 < nz_part[r+1])
                {
                    diag_j[cnt] = row_index + nx_local * ny_local;
                    diag_data[cnt++] = value[6];
                }
                else
                {
                    if (iz+1 < nz)
                    {
                        offd_j[o_cnt] = jxf_map(ix, iy, iz+1, p, q, r+1, P, Q, R,
                                                         nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[6];
                    }
                }
                row_index ++;
            }
        }
    }
    if (num_procs > 1)
    {
        for (i = 0; i < num_cols_offd; i ++)
        {
            col_map_offd[i] = offd_j[i];
        }
        jxf_qsort0(col_map_offd, 0, num_cols_offd-1);
        for (i = 0; i < num_cols_offd; i ++)
        {
            for (j = 0; j < num_cols_offd; j ++)
            {
                if (offd_j[i] == col_map_offd[j])
                {
                    offd_j[i] = j;
                    break;
                }
            }
        }
    }
    A = jxf_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                      global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jxf_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jxf_ParCSRMatrixDiag(A);
    jxf_CSRMatrixI(diag) = diag_i;
    jxf_CSRMatrixJ(diag) = diag_j;
    jxf_CSRMatrixData(diag) = diag_data;
    offd = jxf_ParCSRMatrixOffd(A);
    jxf_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jxf_CSRMatrixJ(offd) = offd_j;
        jxf_CSRMatrixData(offd) = offd_data;
    }
    jxf_TFree(nx_part);
    jxf_TFree(ny_part);
    jxf_TFree(nz_part);
    
    return A;
}
