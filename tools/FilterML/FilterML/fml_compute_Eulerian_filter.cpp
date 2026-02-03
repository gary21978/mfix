#include <AMReX_ParmParse.H>

#include <FilterML.H>
#include <fml_drag_K.H>

using namespace amrex;

void
FilterML::
compute_Eulerian_filter (int const a_filter_size)
{

  Real rho_f, mu_f; // fluid density and viscosity

  { ParmParse pp_fluid("fluid");
    pp_fluid.get("density",  rho_f);
    pp_fluid.get("viscosity", mu_f);
  }

  // number of samples per filter region
  int const nsamples = 3; // af, avg_relV, H_index

  auto drag_coeff = compute_drag_Tang();

  for (int lev(0); lev<get_nlev(); ++lev) {

    const Box& domain = m_geom[lev].Domain();

    BoxArray fml_ba(domain);
    fml_ba.maxSize(a_filter_size);

    BoxList  fml_bl = fml_ba.boxList();

    int const nfilters = static_cast<int>(fml_bl.size());

    int grid_pts, min_gp(INT_MAX), max_gp(0);

    MultiFab const* const p_alpha_f = get_const_alpha_f(lev);
    MultiFab const* const p_alpha_p = get_const_alpha_p(lev);

    MultiFab const* const p_fluid = m_fluid->get_const_data(lev);

    int const Uf_comp = m_fluid->index("Uf");
    int const Vf_comp = m_fluid->index("Vf");
    int const Wf_comp = m_fluid->index("Wf");

    int const dpdy_comp = m_fluid->index("grad_y(Pf)");

    MultiFab const* const p_solids = m_solids->get_const_data(lev);

    int const Up_comp = m_solids->index("Up");
    int const Vp_comp = m_solids->index("Vp");
    int const Wp_comp = m_solids->index("Wp");

    Real const dp = m_particles->get_mean_diameter();

    int const nsums = 15;
    Vector<Real> lev_sums(nsums*nfilters, 0.);

    for (int n(0); n<nfilters; ++n) {

      Box filter_region = fml_ba[n];

      using ROS = amrex::ReduceOpSum;

      ReduceOps<ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS,ROS> reduce_op;
      ReduceData<Real,Real,Real,Real,Real,
                 Real,Real,Real,Real,Real,
                 Real,Real,Real,Real,Real> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

      for (MFIter mfi(*p_alpha_f, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& tbox = mfi.tilebox();

        if ( !tbox.intersects(filter_region)) { continue; }

        Array4<Real const> const& alpha_f  = p_alpha_f->const_array(mfi);

        Array4<Real const> const& Uf  = p_fluid->const_array(mfi, Uf_comp);
        Array4<Real const> const& Vf  = p_fluid->const_array(mfi, Vf_comp);
        Array4<Real const> const& Wf  = p_fluid->const_array(mfi, Wf_comp);
        Array4<Real const> const& dPy = p_fluid->const_array(mfi, dpdy_comp);

        Array4<Real const> const& alpha_p  = p_alpha_p->const_array(mfi);

        Array4<Real const> const& Up = p_solids->const_array(mfi, Up_comp);
        Array4<Real const> const& Vp = p_solids->const_array(mfi, Vp_comp);
        Array4<Real const> const& Wp = p_solids->const_array(mfi, Wp_comp);

        reduce_op.eval(tbox, reduce_data, [filter_region, drag_coeff,
          rho_f, mu_f, alpha_f, Uf, Vf, Wf, dPy, alpha_p, dp, Up, Vp, Wp]
        AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
          Real dcoeff(0), rel_Vel(0.);
          Real af(0.), afUf(0.), afVf(0.), afWf(0.);
          Real ap(0.), apUp(0.), apVp(0.), apWp(0.);

          Real af2(0.), dpdy(0.);
          Real apUf(0.), apVf(0.), apWf(0.);

          if (filter_region.contains(i,j,k)) {

            af = alpha_f(i,j,k);
            afUf = af*Uf(i,j,k);
            afVf = af*Vf(i,j,k);
            afWf = af*Wf(i,j,k);

            ap = alpha_p(i,j,k);
            apUp = ap*Up(i,j,k);
            apVp = ap*Vp(i,j,k);
            apWp = ap*Wp(i,j,k);

            // magnitude of fluid-solids relative velocity
            rel_Vel = std::sqrt(
                (Uf(i,j,k) - Up(i,j,k)) * (Uf(i,j,k) - Up(i,j,k)) +
                (Vf(i,j,k) - Vp(i,j,k)) * (Vf(i,j,k) - Vp(i,j,k)) +
                (Wf(i,j,k) - Wp(i,j,k)) * (Wf(i,j,k) - Wp(i,j,k)));

            dcoeff = (1.-af)*drag_coeff(af, mu_f, rho_f, rel_Vel, dp)*rel_Vel;

            af2 = af*af;
            dpdy = dPy(i,j,k);
            apUf = ap*Uf(i,j,k);
            apVf = ap*Vf(i,j,k);
            apWf = ap*Wf(i,j,k);
          }
          return {af, afUf, afVf, afWf, ap, apUp, apVp, apWp,
                  rel_Vel, dcoeff, af2, dpdy, apUf, apVf, apWf};
        });

      } // MFIter

      ReduceTuple host_tuple = reduce_data.value(reduce_op);

      lev_sums[ 0+n*nsums] = amrex::get<0>(host_tuple);
      lev_sums[ 1+n*nsums] = amrex::get<1>(host_tuple);
      lev_sums[ 2+n*nsums] = amrex::get<2>(host_tuple);
      lev_sums[ 3+n*nsums] = amrex::get<3>(host_tuple);
      lev_sums[ 4+n*nsums] = amrex::get<4>(host_tuple);
      lev_sums[ 5+n*nsums] = amrex::get<5>(host_tuple);
      lev_sums[ 6+n*nsums] = amrex::get<6>(host_tuple);
      lev_sums[ 7+n*nsums] = amrex::get<7>(host_tuple);
      lev_sums[ 8+n*nsums] = amrex::get<8>(host_tuple);
      lev_sums[ 9+n*nsums] = amrex::get<9>(host_tuple);
      lev_sums[10+n*nsums] = amrex::get<10>(host_tuple);
      lev_sums[11+n*nsums] = amrex::get<11>(host_tuple);
      lev_sums[12+n*nsums] = amrex::get<12>(host_tuple);
      lev_sums[13+n*nsums] = amrex::get<12>(host_tuple);
      lev_sums[14+n*nsums] = amrex::get<12>(host_tuple);

      grid_pts = filter_region.numPts();
      min_gp = min(min_gp, grid_pts);
      max_gp = max(max_gp, grid_pts);
    } // loop over filter regions (n)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(min_gp == max_gp,
      "Filter regions not same size!\n");

    ParallelDescriptor::ReduceRealSum(lev_sums.data(), lev_sums.size(),
      ParallelDescriptor::IOProcessorNumber());

    Real const inv_grid_pts = 1.0 / static_cast<Real>(grid_pts);

    Vector<Real> lev_samples(nfilters*nsamples, 0.);

    if(ParallelDescriptor::IOProcessor()) {

    Real const dx(m_geom[lev].CellSize(0));

      std::string filename = csv_filename(lev, a_filter_size);

      std::ofstream File;
      File.open(filename.c_str(), std::ios::out|std::ios::trunc);

      File << "\"x\",\"y\",\"z\",\"avg_af\",\"avg_relV\",\"H_index\",\"var_af\",\"avg_dPdy\",\"drift_U\",\"drift_V\",\"drift_W\"\n";

      for (int n(0); n<nfilters; ++n) {

        Box filter_region = fml_ba[n];
        const int* blo = filter_region.loVect();

        Real const x = dx*static_cast<Real>(blo[0]);
        Real const y = dx*static_cast<Real>(blo[1]);
        Real const z = dx*static_cast<Real>(blo[2]);

        File << std::setprecision(8) << std::scientific << x << ','
             << std::setprecision(8) << std::scientific << y << ','
             << std::setprecision(8) << std::scientific << z << ',';

        Real const af =   lev_sums[0+n*nsums];

        Real const inv_af = (af < 1.e-15) ? 0.0 : 1.0/af;

        // phase averaged fluid velocities:
        Real const Uf = lev_sums[1+n*nsums] * inv_af;
        Real const Vf = lev_sums[2+n*nsums] * inv_af;
        Real const Wf = lev_sums[3+n*nsums] * inv_af;

        Real const ap   = lev_sums[4+n*nsums];

        Real const inv_ap = (ap < 1.e-15) ? 0.0 : 1.0/ap;

        // phase averaged particle velocities:
        Real const Up = lev_sums[5+n*nsums] * inv_ap;
        Real const Vp = lev_sums[6+n*nsums] * inv_ap;
        Real const Wp = lev_sums[7+n*nsums] * inv_ap;

        // fluid volume fraction, relative velocity magnitude
        // and drag coefficient averaged over filter region.
        Real const avg_af = af * inv_grid_pts;
        Real const avg_ap = ap * inv_grid_pts;

        Real const avg_relV   = lev_sums[8+n*nsums] * inv_grid_pts;
        Real const avg_dcoeff = lev_sums[9+n*nsums] * inv_grid_pts;

        // phase averaged fluid-solids relative velocity magnitude
        Real const rel_Vel = std::sqrt(
            (Uf - Up) * (Uf - Up) +
            (Vf - Vp) * (Vf - Vp) +
            (Wf - Wp) * (Wf - Wp));

        // compute drag coefficient over filtered regionusing the average
        // volume fraction and phase averaged relative velocity magnitued
        Real const dcoeff_f =
          (1.-avg_af)*drag_coeff(avg_af, mu_f, rho_f, rel_Vel, dp)*rel_Vel;

        Real const H_index = (std::abs(dcoeff_f) < 1.e-9) ? 0.0 : avg_dcoeff / dcoeff_f;

        Real const var_af   = lev_sums[10+n*nsums]*inv_grid_pts - avg_af*avg_af;

        Real const avg_dPy  = lev_sums[11+n*nsums]*inv_grid_pts;

        Real const drift_U = lev_sums[12+n*nsums]*inv_grid_pts - avg_ap*Uf;
        Real const drift_V = lev_sums[13+n*nsums]*inv_grid_pts - avg_ap*Vf;
        Real const drift_W = lev_sums[14+n*nsums]*inv_grid_pts - avg_ap*Wf;

        File << std::setprecision(8) << std::scientific << avg_af   << ','
             << std::setprecision(8) << std::scientific << avg_relV << ','
             << std::setprecision(8) << std::scientific << H_index  << ','
             << std::setprecision(8) << std::scientific << var_af   << ','
             << std::setprecision(8) << std::scientific << avg_dPy  << ','
             << std::setprecision(8) << std::scientific << drift_U  << ','
             << std::setprecision(8) << std::scientific << drift_V  << ','
             << std::setprecision(8) << std::scientific << drift_W  << '\n';

      }
      File.flush();
      File.close();
    }

    ParallelDescriptor::Barrier();

  }// lev

}
