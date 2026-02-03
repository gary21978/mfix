#ifndef MFIX_DES_K_H_
#define MFIX_DES_K_H_

#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_RealVect.H>
#include <AMReX_Array.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cross_product (const amrex::RealVect& a,
                    const amrex::RealVect& b,
                    amrex::RealVect& c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dot_product (const amrex::RealVect& a,
                         const amrex::RealVect& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interp_level_set (const amrex::RealVect& pos,
                              const int n_refine,
                              amrex::Array4<amrex::Real const> const& phi,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::RealVect scaled_dxi(0.);
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];

    amrex::Real x = (pos[0] - plo[0]) * scaled_dxi[0];
    amrex::Real y = (pos[1] - plo[1]) * scaled_dxi[1];
    amrex::Real z = (pos[2] - plo[2]) * scaled_dxi[2];

    int i = static_cast<int>(amrex::Math::floor(x));
    int j = static_cast<int>(amrex::Math::floor(y));
    int k = static_cast<int>(amrex::Math::floor(z));

    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;

    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    amrex::Real ls_value = 0;
    ls_value += phi(i,   j  , k  ) * wx_lo * wy_lo * wz_lo;
    ls_value += phi(i+1, j  , k  ) * wx_hi * wy_lo * wz_lo;
    ls_value += phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo;
    ls_value += phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo;
    ls_value += phi(i,   j  , k+1) * wx_lo * wy_lo * wz_hi;
    ls_value += phi(i+1, j  , k+1) * wx_hi * wy_lo * wz_hi;
    ls_value += phi(i  , j+1, k+1) * wx_lo * wy_hi * wz_hi;
    ls_value += phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi;

    return ls_value;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void level_set_normal (const amrex::RealVect& pos,
                       const int n_refine,
                       amrex::RealVect& normal,
                       amrex::Array4<amrex::Real const> const& phi,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::RealVect scaled_dxi(0.);
    scaled_dxi[0] = n_refine * dxi[0];
    scaled_dxi[1] = n_refine * dxi[1];
    scaled_dxi[2] = n_refine * dxi[2];

    amrex::Real x = (pos[0] - plo[0]) * scaled_dxi[0];
    amrex::Real y = (pos[1] - plo[1]) * scaled_dxi[1];
    amrex::Real z = (pos[2] - plo[2]) * scaled_dxi[2];

    int i = static_cast<int>(amrex::Math::floor(x));
    int j = static_cast<int>(amrex::Math::floor(y));
    int k = static_cast<int>(amrex::Math::floor(z));

    amrex::Real wx_hi = x - i;
    amrex::Real wy_hi = y - j;
    amrex::Real wz_hi = z - k;

    amrex::Real wx_lo = 1.0 - wx_hi;
    amrex::Real wy_lo = 1.0 - wy_hi;
    amrex::Real wz_lo = 1.0 - wz_hi;

    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    normal[0] -= phi(i,   j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] += phi(i+1, j  , k  ) * scaled_dxi[0] * wy_lo * wz_lo;
    normal[0] -= phi(i,   j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] += phi(i+1, j+1, k  ) * scaled_dxi[0] * wy_hi * wz_lo;
    normal[0] -= phi(i,   j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] += phi(i+1, j  , k+1) * scaled_dxi[0] * wy_lo * wz_hi;
    normal[0] -= phi(i  , j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;
    normal[0] += phi(i+1, j+1, k+1) * scaled_dxi[0] * wy_hi * wz_hi;

    normal[1] -= phi(i,   j  , k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] += phi(i  , j+1, k  ) * scaled_dxi[1] * wx_lo * wz_lo;
    normal[1] -= phi(i+1, j  , k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] += phi(i+1, j+1, k  ) * scaled_dxi[1] * wx_hi * wz_lo;
    normal[1] -= phi(i,   j  , k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] += phi(i  , j+1, k+1) * scaled_dxi[1] * wx_lo * wz_hi;
    normal[1] -= phi(i+1, j  , k+1) * scaled_dxi[1] * wx_hi * wz_hi;
    normal[1] += phi(i+1, j+1, k+1) * scaled_dxi[1] * wx_hi * wz_hi;

    normal[2] -= phi(i  , j  , k  ) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] += phi(i  , j  , k+1) * scaled_dxi[2] * wx_lo * wy_lo;
    normal[2] -= phi(i+1, j  , k  ) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] += phi(i+1, j  , k+1) * scaled_dxi[2] * wx_hi * wy_lo;
    normal[2] -= phi(i,   j+1, k  ) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] += phi(i  , j+1, k+1) * scaled_dxi[2] * wx_lo * wy_hi;
    normal[2] -= phi(i+1, j+1, k  ) * scaled_dxi[2] * wx_hi * wy_hi;
    normal[2] += phi(i+1, j+1, k+1) * scaled_dxi[2] * wx_hi * wy_hi;

    amrex::Real inv_norm = 1.0 / sqrt(dot_product(normal, normal));
    normal[0] *= inv_norm;
    normal[1] *= inv_norm;
    normal[2] *= inv_norm;
}

#endif
