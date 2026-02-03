#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <mfix.H>
#include <mfix_rw.H>
#include <mfix_fluid.H>
#include <mfix_solids.H>
#include <mfix_dem.H>
#include <mfix_pic.H>

#ifdef AMREX_USE_CONDUIT
#include <AMReX_Conduit_Blueprint.H>
#include "conduit_cpp_to_c.hpp"
#endif

#ifdef AMREX_USE_CATALYST
#include "catalyst.hpp"
#endif

#ifdef AMREX_USE_CATALYST
namespace internal {

void EmptyParticleData(const std::vector<std::string>& real_fields,
                       const std::vector<std::string>& int_fields,
                       conduit::Node& node)
{
  const std::string topology_name = "particles";
  // make a dummy node for catalyst
  const std::string& patch_name = amrex::Concatenate("domain_x", 0, 6);
  conduit::Node &patch = node[patch_name];

  patch["state/domain_id"] = 0;

  std::string coordset_name = topology_name + "+coords";
  conduit::Node &n_coords = patch["coordsets"][coordset_name];
  n_coords["type"] = "explicit";

  // create an explicit points topology
  conduit::Node &n_topo = patch["topologies"][topology_name];
  n_topo["coordset"] = coordset_name;
  n_topo["type"] = "unstructured";
  n_topo["elements/shape"] = "point";
  n_topo["elements/connectivity"].set(conduit::DataType::c_int(0));

  n_coords["values/x"].set(conduit::DataType::c_int(0));
  n_coords["values/y"].set(conduit::DataType::c_int(0));
  n_coords["values/z"].set(conduit::DataType::c_int(0));

  conduit::Node &n_fields = patch["fields"];

  // id field
  conduit::Node &n_f_id = n_fields[topology_name + "_id"];
  n_f_id["topology"] = topology_name;
  n_f_id["association"] = "element";
  n_f_id["values"].set(conduit::DataType::c_int(0));

  // cpu field
  conduit::Node &n_f_cpu = n_fields[topology_name + "_cpu"];
  n_f_cpu["topology"] = topology_name;
  n_f_cpu["association"] = "element";
  n_f_cpu["values"].set(conduit::DataType::c_int(0));

  for (auto& rfield : real_fields)
  {
    conduit::Node & n_f = n_fields[rfield];
    n_f["topology"] = topology_name;
    n_f["association"] = "element";
    n_f["values"].set(conduit::DataType::c_int(0));
  }
  for (auto& ifield : int_fields)
  {
    conduit::Node & n_f = n_fields[ifield];
    n_f["topology"] = topology_name;
    n_f["association"] = "element";
    n_f["values"].set(conduit::DataType::c_int(0));
  }
}
} // namespace internal


void MFIXReadWrite::
RunCatalystAdaptor(int nstep, const Real time) const
{
#ifdef AMREX_USE_CONDUIT
  BL_PROFILE("mfix::RunCatalystAdaptor()");

  conduit::Node node;
  auto& state = node["catalyst/state"];
  state["timestep"].set(nstep);
  state["time"].set(time);

  Vector<std::string> fluid_names;
  Vector<std::unique_ptr<MultiFab>> fluid_data(nlev);

  if (fluid.solve()) {

    init_catalyst_fluid(fluid_names, fluid_data);

    auto& channel = node["catalyst/channels/mesh"];
    channel["type"] = "mesh";

#if defined(AMREX_USE_CUDA)
    channel["memoryspace"] = "cuda";
#else
    channel["memoryspace"] = "host";
#endif

    // AMR metadata stays as you have it, under `channel["amr_metadata"]`
    auto& amr = channel["amr_metadata"];
    amr["n_levels"] = nlev;

    // Force data to be an object BEFORE passing to AMReX
    channel["data"].set(conduit::DataType::object());
    auto& meshData = channel["data"];

    Vector<int> level_steps(nlev,nstep);

    amrex::MultiLevelToBlueprint(nlev, amrex::GetVecOfConstPtrs(fluid_data),
        fluid_names, geom, time, level_steps, ref_ratio, meshData);

    // If this is a multi-domain mesh, collapse the first domain
    if (conduit::blueprint::mesh::is_multi_domain(meshData)) {
      // Assume AMReX created domain_0000xx nodes; take the first child
      int ndom = meshData.number_of_children();
      if (ndom > 0) {
        conduit::Node tmp;
        tmp.update(meshData.child(0));  // copy the first domain's mesh

        meshData.reset();      // clear multi-domain container
        meshData.update(tmp);  // now data has coordsets/topologies/fields at its root
      }
    }


  } // end fluid setup


  if ( m_dem.solve() || m_pic.solve() ) {

    Vector<std::string> pc_names_int;
    Vector<std::string> pc_names_real;

    init_catalyst_pc(pc_names_int, pc_names_real);

    auto& channel = node["catalyst/channels/particles"];
    channel["type"] = "mesh";

    auto& particleData = channel["data"];

    conduit::Node amrexParticles;
    amrex::ParticleContainerToBlueprint(*pc,
              pc_names_real, pc_names_int, amrexParticles);

    particleData.update(amrexParticles);

    if(!particleData.dtype().is_object()) {
      internal::EmptyParticleData(pc_names_real, pc_names_int, particleData);
    }

    // collapse multi-domain particle blueprint, same pattern as mesh
    if (conduit::blueprint::mesh::is_multi_domain(particleData)) {
      int ndom = particleData.number_of_children();
      if (ndom > 0) {
        conduit::Node tmp;
        tmp.update(particleData.child(0));  // take the first domain

        particleData.reset();
        particleData.update(tmp);
      }
    }

  } // end particle setup


  // run catalyst
  catalyst_status err = catalyst_execute(conduit::c_node(&node));
  if (err == catalyst_status_ok) {
    Print() << "Catalyst execute was successful\n";
  } else {
    Print() << "Catalyst execute failed\n";
  }

#endif
}


void MFIXReadWrite::
init_catalyst_fluid ( Vector<std::string>& a_names,
                      Vector<std::unique_ptr<MultiFab>>& a_data ) const
{
  BL_PROFILE("MFIXReadWrite::init_catalyst_fluid");

  const int ngrow = 0;
  int ncomp = 5;

  a_names.push_back("u_g");
  a_names.push_back("v_g");
  a_names.push_back("w_g");
  a_names.push_back("ep_g");
  a_names.push_back("volfrac");

  // Temperature in fluid
  if (fluid.solve_enthalpy()) {
    a_names.push_back("T_g");
    ncomp += 1;
  }

  if ( fluid.solve_species()) {
    for (std::string specie: fluid.species_names()) {
      a_names.push_back("Xg_"+specie);
      ncomp += 1;
    }
  }

  AMREX_ALWAYS_ASSERT(a_names.size() == ncomp);

  for (int lev = 0; lev < nlev; ++lev) {

    a_data[lev] = std::make_unique<MultiFab>(grids[lev], dmap[lev],
        ncomp, ngrow,  MFInfo(), *ebfactory[lev]);

    MultiFab::Copy( *(a_data[lev].get()), *(leveldata_const().vel_const(lev)), 0, 0, 1, 0);
    MultiFab::Copy( *(a_data[lev].get()), *(leveldata_const().vel_const(lev)), 1, 1, 1, 0);
    MultiFab::Copy( *(a_data[lev].get()), *(leveldata_const().vel_const(lev)), 2, 2, 1, 0);
    MultiFab::Copy( *(a_data[lev].get()), *(leveldata_const().epf_const(lev)), 0, 3, 1, 0);
    MultiFab::Copy( *(a_data[lev].get()), ebfactory[lev]->getVolFrac(), 0, 4, 1, 0);

    int lc=5;
    if (fluid.solve_enthalpy()) {
      MultiFab::Copy( *(a_data[lev].get()), *(leveldata_const().T_const(lev)), 0, lc, 1, 0);
      lc += 1;
    }

    // Fluid species mass fractions
    if (fluid.solve_species()) {
      MultiFab::Copy( *(a_data[lev].get()), *(leveldata_const().X_const(lev)), 0, lc, fluid.nspecies(), 0);
      lc += fluid.nspecies();
    }

    amrex::EB_set_covered( *(a_data[lev].get()), 0.0);
  }
}


void MFIXReadWrite::
init_catalyst_pc ( Vector<std::string>& a_names_int,
                   Vector<std::string>& a_names_real ) const
{
  BL_PROFILE("MFIXReadWrite::init_catalyst_pc");

  a_names_int.push_back("phase");
  a_names_int.push_back("state");
#if MFIX_POLYDISPERSE
  a_names_int.push_back("ptype");
#endif

  a_names_real.push_back("radius");
  a_names_real.push_back("density");

  a_names_real.push_back("velx");
  a_names_real.push_back("vely");
  a_names_real.push_back("velz");

  if(m_dem.solve()){
    a_names_real.push_back("omegax");
    a_names_real.push_back("omegay");
    a_names_real.push_back("omegaz");
  } else {
    a_names_real.push_back("grad_tau_x");
    a_names_real.push_back("grad_tau_y");
    a_names_real.push_back("grad_tau_z");
  }

  a_names_real.push_back("statwt");
  a_names_real.push_back("dragcoeff");
  a_names_real.push_back("dragx");
  a_names_real.push_back("dragy");
  a_names_real.push_back("dragz");

  a_names_real.push_back("temperature");
}

#endif
