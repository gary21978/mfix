#ifndef AMREX_SPMV_H_
#define AMREX_SPMV_H_
#include <AMReX_Config.H>

#include <AMReX_AlgVector.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_SpMatrix.H>

#if defined(AMREX_USE_CUDA)
#  include <cusparse.h>
#endif

namespace amrex {

template <typename T>
void SpMV (AlgVector<T>& y, SpMatrix<T> const& A, AlgVector<T> const& x)
{
    // xxxxx TODOL We might want to cache the cusparse and rocsparse handles

    // xxxxx TODO: let's assume it's square matrix for now.
    AMREX_ALWAYS_ASSERT(x.partition() == y.partition() &&
                        x.partition() == A.partition());

    const_cast<SpMatrix<T>&>(A).startComm(x);

    T      * AMREX_RESTRICT py = y.data();
    T const* AMREX_RESTRICT px = x.data();
    T const* AMREX_RESTRICT mat = A.data();
    auto const* AMREX_RESTRICT col = A.columnIndex();
    auto const* AMREX_RESTRICT row = A.rowOffset();

#if defined(AMREX_USE_GPU)

    Long const nrows = A.numLocalRows();
    Long const ncols = x.numLocalRows();
    Long const nnz = A.numLocalNonZero();

#if defined(AMREX_USE_CUDA)

    cusparseHandle_t handle;
    cusparseCreate(&handle);
    cusparseSetStream(handle, Gpu::gpuStream());

    cudaDataType data_type;
    if constexpr (std::is_same_v<T,float>) {
        data_type = CUDA_R_32F;
    } else if constexpr (std::is_same_v<T,double>) {
        data_type = CUDA_R_64F;
    } else if constexpr (std::is_same_v<T,GpuComplex<float>>) {
        data_type = CUDA_C_32F;
    } else if constexpr (std::is_same_v<T,GpuComplex<double>>) {
        data_type = CUDA_C_64F;
    } else {
        amrex::Abort("SpMV: unsupported data type");
    }

    cusparseIndexType_t index_type = CUSPARSE_INDEX_64I;

    cusparseSpMatDescr_t mat_descr;
    cusparseCreateCsr(&mat_descr, nrows, ncols, nnz, (void*)row, (void*)col, (void*)mat,
                      index_type, index_type, CUSPARSE_INDEX_BASE_ZERO, data_type);

    cusparseDnVecDescr_t x_descr;
    cusparseCreateDnVec(&x_descr, ncols, (void*)px, data_type);

    cusparseDnVecDescr_t y_descr;
    cusparseCreateDnVec(&y_descr, nrows, (void*)py, data_type);

    T alpha = T(1);
    T beta = T(0);

    std::size_t buffer_size;
    cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_descr, x_descr,
                            &beta, y_descr, data_type, CUSPARSE_SPMV_ALG_DEFAULT, &buffer_size);

    auto* pbuffer = (void*)The_Arena()->alloc(buffer_size);

    cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_descr, x_descr,
                 &beta, y_descr, data_type, CUSPARSE_SPMV_ALG_DEFAULT, pbuffer);

    Gpu::streamSynchronize();

    cusparseDestroySpMat(mat_descr);
    cusparseDestroyDnVec(x_descr);
    cusparseDestroyDnVec(y_descr);
    cusparseDestroy(handle);
    The_Arena()->free(pbuffer);

#else
    amrex::Abort("SpMV: CUDA backend required for GPU path");
#endif

    AMREX_GPU_ERROR_CHECK();

#else

    Long const ny = y.numLocalRows();

    for (Long i = 0; i < ny; ++i) {
        T r = 0;
        for (Long j = row[i]; j < row[i+1]; ++j) {
            r += mat[j] * px[col[j]];
        }
        py[i] = r;
    }

#endif

    const_cast<SpMatrix<T>&>(A).finishComm(y);
}

}

#endif
