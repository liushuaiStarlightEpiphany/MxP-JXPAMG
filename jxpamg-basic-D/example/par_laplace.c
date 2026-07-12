//========================================================================//
//  JXPAMG(IAPCM & XTU Parallel Algebraic Multigrid) (c) 2009-2013        //
//  Institute of Applied Physics and Computational Mathematics            //
//  School of Mathematics and Computational Science Xiangtan University   //
//========================================================================//

/*!
 *  par_laplace.c -- Laplace problem
 *  Date: 2012/09/19
 *
 *  Created by Yue Xiaoqiang
 */

#include "jx_util.h"
#include "jx_mv.h"

/*!
 * \fn JX_Int jx_map
 */
JX_Int
jx_map( JX_Int ix,
        JX_Int iy,
        JX_Int iz,
        JX_Int p,
        JX_Int q,
        JX_Int r,
        JX_Int P,
        JX_Int Q,
        JX_Int R,
        JX_Int *nx_part,
        JX_Int *ny_part,
        JX_Int *nz_part,
        JX_Int *global_part )
{
    JX_Int nx_local;
    JX_Int ny_local;
    JX_Int ix_local;
    JX_Int iy_local;
    JX_Int iz_local;
    JX_Int global_index;
    JX_Int proc_num;
    
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
 * \fn JX_Int jx_map2
 */
JX_Int
jx_map2( JX_Int ix,
         JX_Int iy,
         JX_Int p,
         JX_Int q,
         JX_Int P,
         JX_Int Q,
         JX_Int *nx_part,
         JX_Int *ny_part,
         JX_Int *global_part )
{
    JX_Int nx_local;
    JX_Int ix_local;
    JX_Int iy_local;
    JX_Int global_index;
    JX_Int proc_num;
    
    proc_num = q * P + p;
    nx_local = nx_part[p+1] - nx_part[p];
    ix_local = ix - nx_part[p];
    iy_local = iy - ny_part[q];
    global_index = global_part[proc_num] + iy_local * nx_local + ix_local;
    
    return global_index;
}

/*!
 * \fn JX_Int jx_map3
 */
JX_Int
jx_map3( JX_Int ix,
         JX_Int iy,
         JX_Int iz,
         JX_Int p,
         JX_Int q,
         JX_Int r,
         JX_Int P,
         JX_Int Q,
         JX_Int R,
         JX_Int *nx_part,
         JX_Int *ny_part,
         JX_Int *nz_part,
         JX_Int *global_part )
{
    JX_Int nx_local;
    JX_Int ix_local;
    JX_Int iy_local;
    JX_Int iz_local;
    JX_Int nxy;
    JX_Int global_index;
    JX_Int proc_num;
    
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
 * \fn jx_ParCSRMatrix *jx_GenerateParLaplacian
 */
jx_ParCSRMatrix *
jx_GenerateParLaplacian( MPI_Comm comm,
                         JX_Int nx,
                         JX_Int ny,
                         JX_Int nz,
                         JX_Int P,
                         JX_Int Q,
                         JX_Int R,
                         JX_Int p,
                         JX_Int q,
                         JX_Int r,
                         JX_Real *value )
{
    jx_ParCSRMatrix *A = NULL;
    
    jx_CSRMatrix *diag = NULL;
    jx_CSRMatrix *offd = NULL;
    
    JX_Int *diag_i = NULL;
    JX_Int *diag_j = NULL;
    JX_Real *diag_data = NULL;
    
    JX_Int *offd_i = NULL;
    JX_Int *offd_j = NULL;
    JX_Real *offd_data = NULL;
    
    JX_Int *global_part = NULL;
    JX_Int *col_map_offd = NULL;
    
    JX_Int *nx_part = NULL;
    JX_Int *ny_part = NULL;
    JX_Int *nz_part = NULL;
    
    JX_Int ix, iy, iz;
    JX_Int cnt, o_cnt;
    JX_Int local_num_rows;
    JX_Int row_index;
    JX_Int i, j;
    
    JX_Int nx_local, ny_local, nz_local;
    JX_Int nx_size, ny_size, nz_size;
    JX_Int num_cols_offd;
    JX_Int grid_size;
    
    JX_Int num_procs, my_id;
    JX_Int P_busy, Q_busy, R_busy;
    
    jx_MPI_Comm_size(comm, &num_procs);
    jx_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny * nz;
    
    jx_GeneratePartitioning(nx, P, &nx_part);
    jx_GeneratePartitioning(ny, Q, &ny_part);
    jx_GeneratePartitioning(nz, R, &nz_part);
    
    global_part = jx_CTAlloc(JX_Int, P*Q*R+1);
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
    
    diag_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    offd_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    P_busy = jx_min(nx, P);
    Q_busy = jx_min(ny, Q);
    R_busy = jx_min(nz, R);
    
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
    
    col_map_offd = jx_CTAlloc(JX_Int, num_cols_offd);
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
    diag_j = jx_CTAlloc(JX_Int, diag_i[local_num_rows]);
    diag_data = jx_CTAlloc(JX_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jx_CTAlloc(JX_Int, offd_i[local_num_rows]);
        offd_data = jx_CTAlloc(JX_Real, offd_i[local_num_rows]);
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
                        offd_j[o_cnt] = jx_map(ix, iy, iz-1, p, q, r-1, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix, iy-1, iz, p, q-1, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix-1, iy, iz, p-1, q, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix+1, iy, iz, p+1, q, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix, iy+1, iz, p, q+1, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix, iy, iz+1, p, q, r+1, P, Q, R,
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
        jx_qsort0(col_map_offd, 0, num_cols_offd-1);
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
#ifdef JX_NO_GLOBAL_PARTITION
/* ideally we would use less storage earlier in this function, but this is fine for testing */
    {
        JX_Int tmp1, tmp2;
        tmp1 = global_part[my_id];
        tmp2 = global_part[my_id+1];
        jx_TFree(global_part);
        global_part = jx_CTAlloc(JX_Int, 2);
        global_part[0] = tmp1;
        global_part[1] = tmp2;
    }
#endif
    A = jx_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                      global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jx_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jx_ParCSRMatrixDiag(A);
    jx_CSRMatrixI(diag) = diag_i;
    jx_CSRMatrixJ(diag) = diag_j;
    jx_CSRMatrixData(diag) = diag_data;
    offd = jx_ParCSRMatrixOffd(A);
    jx_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jx_CSRMatrixJ(offd) = offd_j;
        jx_CSRMatrixData(offd) = offd_data;
    }
    jx_TFree(nx_part);
    jx_TFree(ny_part);
    jx_TFree(nz_part);
    
    return A;
}

/*!
 * \fn jx_ParCSRMatrix *jx_GenerateParLaplacian2d9pt
 */
jx_ParCSRMatrix *
jx_GenerateParLaplacian2d9pt( MPI_Comm comm,
                              JX_Int nx,
                              JX_Int ny,
                              JX_Int P,
                              JX_Int Q,
                              JX_Int p,
                              JX_Int q,
                              JX_Real *value )
{
    jx_ParCSRMatrix *A = NULL;
    
    jx_CSRMatrix *diag = NULL;
    jx_CSRMatrix *offd = NULL;
    
    JX_Int *diag_i = NULL;
    JX_Int *diag_j = NULL;
    JX_Real *diag_data = NULL;
    
    JX_Int *offd_i = NULL;
    JX_Int *offd_j = NULL;
    JX_Real *offd_data = NULL;
    
    JX_Int *work = NULL;
    JX_Int *nx_part = NULL;
    JX_Int *ny_part = NULL;
    JX_Int *global_part = NULL;
    JX_Int *col_map_offd = NULL;
    
    JX_Int i, j, ix, iy, cnt, o_cnt;
    JX_Int num_cols_offd, grid_size;
    JX_Int local_num_rows, row_index;
    JX_Int num_procs, my_id, P_busy, Q_busy;
    JX_Int nx_local, ny_local, nx_size, ny_size;
    
    jx_MPI_Comm_size(comm, &num_procs);
    jx_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny;
    
    jx_GeneratePartitioning(nx, P, &nx_part);
    jx_GeneratePartitioning(ny, Q, &ny_part);
    
    global_part = jx_CTAlloc(JX_Int, P*Q+1);
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
    
    diag_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    offd_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    P_busy = jx_min(nx, P);
    Q_busy = jx_min(ny, Q);
    
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
    col_map_offd = jx_CTAlloc(JX_Int, num_cols_offd);
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
    diag_j = jx_CTAlloc(JX_Int, diag_i[local_num_rows]);
    diag_data = jx_CTAlloc(JX_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jx_CTAlloc(JX_Int, offd_i[local_num_rows]);
        offd_data = jx_CTAlloc(JX_Real, offd_i[local_num_rows]);
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
                        offd_j[o_cnt] = jx_map2(ix-1, iy-1, p-1, q, P, Q, nx_part, ny_part, global_part);
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
                        offd_j[o_cnt] = jx_map2(ix+1, iy-1, p+1, q, P, Q, nx_part, ny_part, global_part);
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
                        offd_j[o_cnt] = jx_map2(ix-1, iy-1, p, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix)
                    {
                        offd_j[o_cnt] = jx_map2(ix-1, iy-1, p-1, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    offd_j[o_cnt] = jx_map2(ix, iy-1, p, q-1, P, Q, nx_part, ny_part, global_part);
                    offd_data[o_cnt++] = value[1];
                    if (ix < nx_part[p+1]-1)
                    {
                        offd_j[o_cnt] = jx_map2(ix+1, iy-1, p, q-1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix+1 < nx)
                    {
                        offd_j[o_cnt] = jx_map2(ix+1, iy-1, p+1, q-1, P, Q, nx_part, ny_part, global_part);
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
                    offd_j[o_cnt] = jx_map2(ix-1, iy, p-1, q, P, Q, nx_part, ny_part, global_part);
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
                    offd_j[o_cnt] = jx_map2(ix+1, iy, p+1, q, P, Q, nx_part, ny_part, global_part);
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
                        offd_j[o_cnt] = jx_map2(ix-1, iy+1, p-1, q, P, Q, nx_part, ny_part, global_part);
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
                        offd_j[o_cnt] = jx_map2(ix+1, iy+1, p+1, q, P, Q, nx_part, ny_part, global_part);
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
                        offd_j[o_cnt] = jx_map2(ix-1, iy+1, p, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix)
                    {
                        offd_j[o_cnt] = jx_map2(ix-1, iy+1, p-1, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    offd_j[o_cnt] = jx_map2(ix, iy+1, p, q+1, P, Q, nx_part, ny_part, global_part);
                    offd_data[o_cnt++] = value[1];
                    if (ix < nx_part[p+1]-1)
                    {
                        offd_j[o_cnt] = jx_map2(ix+1, iy+1, p, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                    else if (ix < nx-1)
                    {
                        offd_j[o_cnt] = jx_map2(ix+1, iy+1, p+1, q+1, P, Q, nx_part, ny_part, global_part);
                        offd_data[o_cnt++] = value[1];
                    }
                }
            }
            row_index ++;
        }
    }
    if (num_procs > 1)
    {
        work = jx_CTAlloc(JX_Int, o_cnt);
        for (i = 0; i < o_cnt; i ++)
        {
            work[i] = offd_j[i];
        }
        jx_qsort0(work, 0, o_cnt-1);
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
        jx_TFree(work);
    }
#ifdef JX_NO_GLOBAL_PARTITION
/* ideally we would use less storage earlier in this function, but this is fine for testing */
    {
        JX_Int tmp1, tmp2;
        tmp1 = global_part[my_id];
        tmp2 = global_part[my_id+1];
        jx_TFree(global_part);
        global_part = jx_CTAlloc(JX_Int, 2);
        global_part[0] = tmp1;
        global_part[1] = tmp2;
    }
#endif
    A = jx_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                 global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jx_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jx_ParCSRMatrixDiag(A);
    jx_CSRMatrixI(diag) = diag_i;
    jx_CSRMatrixJ(diag) = diag_j;
    jx_CSRMatrixData(diag) = diag_data;
    offd = jx_ParCSRMatrixOffd(A);
    jx_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jx_CSRMatrixJ(offd) = offd_j;
        jx_CSRMatrixData(offd) = offd_data;
    }
    jx_TFree(nx_part);
    jx_TFree(ny_part);
    
    return A;
}

/*!
 * \fn jx_ParCSRMatrix *jx_GenerateParLaplacian2d9pt
 */
jx_ParCSRMatrix *
jx_GenerateParLaplacian3d27pt( MPI_Comm comm,
                               JX_Int nx,
                               JX_Int ny,
                               JX_Int nz,
                               JX_Int P,
                               JX_Int Q,
                               JX_Int R,
                               JX_Int p,
                               JX_Int q,
                               JX_Int r,
                               JX_Real *value )
{
    jx_ParCSRMatrix *A = NULL;
    
    jx_CSRMatrix *diag = NULL;
    jx_CSRMatrix *offd = NULL;
    
    JX_Int *diag_i = NULL;
    JX_Int *diag_j = NULL;
    JX_Real *diag_data = NULL;
    
    JX_Int *offd_i = NULL;
    JX_Int *offd_j = NULL;
    JX_Real *offd_data = NULL;
    
    JX_Int *work = NULL;
    JX_Int *nx_part = NULL;
    JX_Int *ny_part = NULL;
    JX_Int *nz_part = NULL;
    JX_Int *global_part = NULL;
    JX_Int *col_map_offd = NULL;
    
    JX_Int i, j, ix, iy, iz;
    JX_Int num_cols_offd, nxy, grid_size;
    JX_Int cnt, o_cnt, local_num_rows, row_index;
    JX_Int num_procs, my_id, P_busy, Q_busy, R_busy;
    JX_Int nx_local, ny_local, nz_local, nx_size, ny_size, nz_size;
    
    jx_MPI_Comm_size(comm, &num_procs);
    jx_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny * nz;
    
    jx_GeneratePartitioning(nx, P, &nx_part);
    jx_GeneratePartitioning(ny, Q, &ny_part);
    jx_GeneratePartitioning(nz, R, &nz_part);
    
    global_part = jx_CTAlloc(JX_Int, P*Q*R+1);
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
    diag_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    offd_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    
    P_busy = jx_min(nx, P);
    Q_busy = jx_min(ny, Q);
    R_busy = jx_min(nz, R);
    
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
    col_map_offd = jx_CTAlloc(JX_Int, num_cols_offd);
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
    diag_j = jx_CTAlloc(JX_Int, diag_i[local_num_rows]);
    diag_data = jx_CTAlloc(JX_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jx_CTAlloc(JX_Int, offd_i[local_num_rows]);
        offd_data = jx_CTAlloc(JX_Real, offd_i[local_num_rows]);
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p-1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p+1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p-1, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy-1, iz-1, p, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p+1, q-1, r,
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
                            offd_j[o_cnt] = jx_map3(ix-1, iy, iz-1, p-1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix+1, iy, iz-1, p+1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p-1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p+1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p-1, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy+1, iz-1, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p+1, q+1, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p-1, q, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy-1, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p+1, q, r-1,
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
                                    offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz-1, p-1, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jx_map3(ix, iy-1, iz-1, p, q-1, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz-1, p+1, q-1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        if (ix > nx_part[p])
                        {
                            offd_j[o_cnt] = jx_map3(ix-1, iy, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy, iz-1, p-1, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        offd_j[o_cnt] = jx_map3(ix, iy, iz-1, p, q, r-1,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix+1 < nx_part[p+1])
                        {
                            offd_j[o_cnt] = jx_map3(ix+1, iy, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy, iz-1, p+1, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        if (iy+1 < ny_part[q+1])
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p-1, q, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy+1, iz-1, p, q, r-1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p, q, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p+1, q, r-1,
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
                                    offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz-1, p-1, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jx_map3(ix, iy+1, iz-1, p, q+1, r-1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p, q+1, r-1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz-1, p+1, q+1, r-1,
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
                            offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz, p-1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz, p+1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz, p, q-1, r,
                                               P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix)
                        {
                            offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz, p-1, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        offd_j[o_cnt] = jx_map3(ix, iy-1, iz, p, q-1, r,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz, p, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix < nx-1)
                        {
                            offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz, p+1, q-1, r,
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
                        offd_j[o_cnt] = jx_map3(ix-1, iy, iz, p-1, q, r,
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
                        offd_j[o_cnt] = jx_map3(ix+1, iy, iz, p+1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz, p-1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz, p+1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix)
                        {
                            offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz, p-1, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        offd_j[o_cnt] = jx_map3(ix, iy+1, iz, p, q+1, r,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix < nx_part[p+1]-1)
                        {
                            offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else if (ix < nx-1)
                        {
                            offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz, p+1, q+1, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p-1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p+1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p-1, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy-1, iz+1, p, q-1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p, q-1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p+1, q-1, r,
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
                            offd_j[o_cnt] = jx_map3(ix-1, iy, iz+1, p-1, q, r,
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
                            offd_j[o_cnt] = jx_map3(ix+1, iy, iz+1, p+1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p-1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p+1, q, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix)
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p-1, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy+1, iz+1, p, q+1, r,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p, q+1, r,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else if (ix < nx-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p+1, q+1, r,
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
                                offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p-1, q, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy-1, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p+1, q, r+1,
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
                                    offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy-1, iz+1, p-1, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jx_map3(ix, iy-1, iz+1, p, q-1, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy-1, iz+1, p+1, q-1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                        }
                        if (ix > nx_part[p])
                        {
                            offd_j[o_cnt] = jx_map3(ix-1, iy, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix)
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy, iz+1, p-1, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        offd_j[o_cnt] = jx_map3(ix, iy, iz+1, p, q, r+1,
                                            P, Q, R, nx_part, ny_part, nz_part, global_part);
                        offd_data[o_cnt++] = value[1];
                        if (ix+1 < nx_part[p+1])
                        {
                            offd_j[o_cnt] = jx_map3(ix+1, iy, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                        }
                        else
                        {
                            if (ix+1 < nx)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy, iz+1, p+1, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                        }
                        if (iy+1 < ny_part[q+1])
                        {
                            if (ix > nx_part[p])
                            {
                                offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p-1, q, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                            }
                            offd_j[o_cnt] = jx_map3(ix, iy+1, iz+1, p, q, r+1,
                                                P, Q, R, nx_part, ny_part, nz_part, global_part);
                            offd_data[o_cnt++] = value[1];
                            if (ix < nx_part[p+1]-1)
                            {
                                offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p, q, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                            }
                            else
                            {
                                if (ix+1 < nx)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p+1, q, r+1,
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
                                    offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix)
                                {
                                    offd_j[o_cnt] = jx_map3(ix-1, iy+1, iz+1, p-1, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                offd_j[o_cnt] = jx_map3(ix, iy+1, iz+1, p, q+1, r+1,
                                                    P, Q, R, nx_part, ny_part, nz_part, global_part);
                                offd_data[o_cnt++] = value[1];
                                if (ix < nx_part[p+1]-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p, q+1, r+1,
                                                        P, Q, R, nx_part, ny_part, nz_part, global_part);
                                    offd_data[o_cnt++] = value[1];
                                }
                                else if (ix < nx-1)
                                {
                                    offd_j[o_cnt] = jx_map3(ix+1, iy+1, iz+1, p+1, q+1, r+1,
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
        work = jx_CTAlloc(JX_Int, o_cnt);
        for (i = 0; i < o_cnt; i ++)
        {
            work[i] = offd_j[i];
        }
        jx_qsort0(work, 0, o_cnt-1);
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
        jx_TFree(work);
    }
#ifdef JX_NO_GLOBAL_PARTITION
/* ideally we would use less storage earlier in this function, but this is fine for testing */
    {
        JX_Int tmp1, tmp2;
        tmp1 = global_part[my_id];
        tmp2 = global_part[my_id+1];
        jx_TFree(global_part);
        global_part = jx_CTAlloc(JX_Int, 2);
        global_part[0] = tmp1;
        global_part[1] = tmp2;
    }
#endif
    A = jx_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                   global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jx_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jx_ParCSRMatrixDiag(A);
    jx_CSRMatrixI(diag) = diag_i;
    jx_CSRMatrixJ(diag) = diag_j;
    jx_CSRMatrixData(diag) = diag_data;
    offd = jx_ParCSRMatrixOffd(A);
    jx_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jx_CSRMatrixJ(offd) = offd_j;
        jx_CSRMatrixData(offd) = offd_data;
    }
    jx_TFree(nx_part);
    jx_TFree(ny_part);
    jx_TFree(nz_part);
    
    return A;
}

/*!
 * \fn jx_ParCSRMatrix *jx_GenerateParConvecDiff
 */
jx_ParCSRMatrix *
jx_GenerateParConvecDiff( MPI_Comm comm,
                          JX_Int nx,
                          JX_Int ny,
                          JX_Int nz,
                          JX_Int P,
                          JX_Int Q,
                          JX_Int R,
                          JX_Int p,
                          JX_Int q,
                          JX_Int r,
                          JX_Real  *value )
{
    jx_ParCSRMatrix *A = NULL;
    
    jx_CSRMatrix *diag = NULL;
    jx_CSRMatrix *offd = NULL;
    
    JX_Int *diag_i = NULL;
    JX_Int *diag_j = NULL;
    JX_Real *diag_data = NULL;
    
    JX_Int *offd_i = NULL;
    JX_Int *offd_j = NULL;
    JX_Real *offd_data = NULL;
    
    JX_Int *nx_part = NULL;
    JX_Int *ny_part = NULL;
    JX_Int *nz_part = NULL;
    JX_Int *global_part = NULL;
    JX_Int *col_map_offd = NULL;
    
    JX_Int i, j, ix, iy, iz;
    JX_Int num_cols_offd, grid_size;
    JX_Int row_index, cnt, o_cnt, local_num_rows;
    JX_Int num_procs, my_id, P_busy, Q_busy, R_busy;
    JX_Int nx_local, ny_local, nz_local, nx_size, ny_size, nz_size;
    
    jx_MPI_Comm_size(comm, &num_procs);
    jx_MPI_Comm_rank(comm, &my_id);
    
    grid_size = nx * ny * nz;
    
    jx_GeneratePartitioning(nx, P, &nx_part);
    jx_GeneratePartitioning(ny, Q, &ny_part);
    jx_GeneratePartitioning(nz, R, &nz_part);
    
    global_part = jx_CTAlloc(JX_Int, P*Q*R+1);
    
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
    diag_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    offd_i = jx_CTAlloc(JX_Int, local_num_rows+1);
    
    P_busy = jx_min(nx, P);
    Q_busy = jx_min(ny, Q);
    R_busy = jx_min(nz, R);
    
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
    col_map_offd = jx_CTAlloc(JX_Int, num_cols_offd);
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
    diag_j = jx_CTAlloc(JX_Int, diag_i[local_num_rows]);
    diag_data = jx_CTAlloc(JX_Real, diag_i[local_num_rows]);
    if (num_procs > 1)
    {
        offd_j = jx_CTAlloc(JX_Int, offd_i[local_num_rows]);
        offd_data = jx_CTAlloc(JX_Real, offd_i[local_num_rows]);
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
                        offd_j[o_cnt] = jx_map(ix, iy, iz-1, p, q, r-1, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix, iy-1, iz, p, q-1, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix-1, iy, iz, p-1, q, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix+1, iy, iz, p+1, q, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix, iy+1, iz, p, q+1, r, P, Q, R,
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
                        offd_j[o_cnt] = jx_map(ix, iy, iz+1, p, q, r+1, P, Q, R,
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
        jx_qsort0(col_map_offd, 0, num_cols_offd-1);
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
    A = jx_ParCSRMatrixCreate(comm, grid_size, grid_size, global_part,
                      global_part, num_cols_offd, diag_i[local_num_rows], offd_i[local_num_rows]);
    jx_ParCSRMatrixColMapOffd(A) = col_map_offd;
    diag = jx_ParCSRMatrixDiag(A);
    jx_CSRMatrixI(diag) = diag_i;
    jx_CSRMatrixJ(diag) = diag_j;
    jx_CSRMatrixData(diag) = diag_data;
    offd = jx_ParCSRMatrixOffd(A);
    jx_CSRMatrixI(offd) = offd_i;
    if (num_cols_offd)
    {
        jx_CSRMatrixJ(offd) = offd_j;
        jx_CSRMatrixData(offd) = offd_data;
    }
    jx_TFree(nx_part);
    jx_TFree(ny_part);
    jx_TFree(nz_part);
    
    return A;
}
