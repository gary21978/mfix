#include <mfix.H>

#include <mfix_mf_helpers.H>
#include <mfix_eb.H>
#include <mfix_dem.H>
#include <mfix_fluid.H>
#include <mfix_species.H>
#include <mfix_reactions.H>
#include <mfix_pic.H>
#include <mfix_monitors.H>

#include <mfix_stepdata.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

void mfix::
EvolveFluid ( MFIXStepData& /*a_stepData*/,
              int /*nstep*/,
              Real& /*dt*/,
              Real& /*prev_dt*/,
              const Real /*time*/,
              Real /*stop_time*/,
              Real& /*coupling_timing*/ )
{
// Saving the empty file because some of the functionality
// currently in the evolve routine may get relocated
// back here once the reorganization is finished.
}
