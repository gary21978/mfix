#ifndef _MFIX_EB_NEWTON_SOLVER_K_H_
#define _MFIX_EB_NEWTON_SOLVER_K_H_

/****************************************************************
 *                                                              *
 *                                                              *
 *                                                              *
 *                                                              *
 ***************************************************************/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real f (const int dir,
               const amrex::RealVect& pos,
               const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
               const amrex::Real& xi,
               const amrex::Real& eta,
               const amrex::Real& zeta)
{
    amrex::Real a0 = (nodes[0][dir] - pos[dir]);
    amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
    amrex::Real a2 = (nodes[2][dir] - nodes[0][dir]);
    amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a0 + a1*xi + a2*eta + a3*zeta + a4*xi*eta + a5*xi*zeta + a6*eta*zeta + a7*xi*eta*zeta;
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dfdxi (const int dir,
                   const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                   const amrex::Real& /*xi*/,
                   const amrex::Real& eta,
                   const amrex::Real& zeta)
{
    amrex::Real a1 = (nodes[1][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a1 + a4*eta + a5*zeta + a7*eta*zeta;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dfdeta (const int dir,
                    const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                    const amrex::Real& xi,
                    const amrex::Real& /*eta*/,
                    const amrex::Real& zeta)
{

    amrex::Real a2 = (nodes[2][dir] - nodes[0][dir]);
    amrex::Real a4 = (nodes[0][dir] - nodes[1][dir] + nodes[3][dir] - nodes[2][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a2 + a4*xi + a6*zeta + a7*xi*zeta;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real dfdzeta (const int dir,
                     const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                     const amrex::Real& xi,
                     const amrex::Real& eta,
                     const amrex::Real& /*zeta*/)
{

    amrex::Real a3 = (nodes[4][dir] - nodes[0][dir]);
    amrex::Real a5 = (nodes[0][dir] - nodes[1][dir] - nodes[4][dir] + nodes[5][dir]);
    amrex::Real a6 = (nodes[0][dir] - nodes[2][dir] - nodes[4][dir] + nodes[6][dir]);
    amrex::Real a7 = (nodes[1][dir] - nodes[3][dir] + nodes[2][dir] + nodes[4][dir] - nodes[5][dir] + nodes[7][dir] - nodes[6][dir] - nodes[0][dir]);

    return a3 + a5*xi + a6*eta + a7*xi*eta;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_interp_mapping (const amrex::RealVect& pos,
                         const amrex::GpuArray<amrex::GpuArray<amrex::Real,3>,8>& nodes,
                         amrex::Real& xi,
                         amrex::Real& eta,
                         amrex::Real& zeta)
{


    int lc(0);
    amrex::Real err(1.0);

    while(err > 1.0e-3 && lc < 10) {

        amrex::Real f0 = f(0, pos, nodes, xi, eta, zeta);
        amrex::Real f1 = f(1, pos, nodes, xi, eta, zeta);
        amrex::Real f2 = f(2, pos, nodes, xi, eta, zeta);

        amrex::Real df0dxi   = dfdxi(0, nodes, xi, eta, zeta);
        amrex::Real df0deta  = dfdeta(0, nodes, xi, eta, zeta);
        amrex::Real df0dzeta = dfdzeta(0, nodes, xi, eta, zeta);

        amrex::Real df1dxi   = dfdxi(1, nodes, xi, eta, zeta);
        amrex::Real df1deta  = dfdeta(1, nodes, xi, eta, zeta);
        amrex::Real df1dzeta = dfdzeta(1, nodes, xi, eta, zeta);

        amrex::Real df2dxi   = dfdxi(2, nodes, xi, eta, zeta);
        amrex::Real df2deta  = dfdeta(2, nodes, xi, eta, zeta);
        amrex::Real df2dzeta = dfdzeta(2, nodes, xi, eta, zeta);

        amrex::Real detJ = df0dxi   * (df1deta*df2dzeta - df1dzeta*df2deta)
                         - df0deta  * (df1dxi*df2dzeta  - df1dzeta*df2dxi)
                         + df0dzeta * (df1dxi*df2deta   - df1deta*df2dxi);

        amrex::Real detJ_xi = f0       * (df1deta*df2dzeta - df1dzeta*df2deta)
                            - df0deta  * (f1*df2dzeta      - df1dzeta*f2)
                            + df0dzeta * (f1*df2deta       - df1deta*f2);

        amrex::Real detJ_eta = df0dxi   * (f1*df2dzeta     - df1dzeta*f2)
                             - f0       * (df1dxi*df2dzeta - df1dzeta*df2dxi)
                             + df0dzeta * (df1dxi*f2       - f1*df2dxi);

        amrex::Real detJ_zeta = df0dxi   * (df1deta*f2     - f1*df2deta)
                              - df0deta  * (df1dxi*f2      - f1*df2dxi)
                              + f0       * (df1dxi*df2deta - df1deta*df2dxi);

        amrex::Real new_xi   = xi   - detJ_xi   / detJ;
        amrex::Real new_eta  = eta  - detJ_eta  / detJ;
        amrex::Real new_zeta = zeta - detJ_zeta / detJ;

        lc += 1;
        err = amrex::max(amrex::Math::abs(xi-new_xi), amrex::Math::abs(eta-new_eta), amrex::Math::abs(zeta-new_zeta));

        xi   = new_xi;
        eta  = new_eta;
        zeta = new_zeta;

    }

}

#endif
