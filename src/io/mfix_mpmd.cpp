#include <AMReX_Print.H>
#include <mfix_rw.H>
#include <sstream>

using namespace amrex;

void MFIXReadWrite::ConstructMPMDHeader (std::string& header) const
{
  std::stringstream ss;

  ss << "{\"data\": {\"static_mfs\": [";

  for (int i = 0; i < m_mpmd_static_mf_names.size(); ++i) {
    int c = 1;
    if (m_mpmd_static_mf_names[i] == "centroid") {
      c = 3;
    }
    ss << "{\"n\": \"" << m_mpmd_static_mf_names[i] << "\", \"c\":" << c << "}";
    if (i < (m_mpmd_static_mf_names.size() - 1)) {
      ss << ",";
    }
  }

  ss << "], \"int_flags_root\": [\"end\"], \"reals_root\": [\"time\"], \"mfs\": [";

  Vector<std::string> names;
  Vector<int> comps;
  for (const auto& name: m_mpmd_mf_names) {
    if ( (name == "ep_g") || (name == "T_g") ) {
      names.push_back(name);
      comps.push_back(1);
    }
    else if (name == "vel_g") {
      names.push_back(name);
      comps.push_back(3);
    } else if (name == "X_gk") {
      for (const auto& specie: fluid.species_names()) {
        names.push_back("X_"+specie+"_g");
        comps.push_back(1);
      }
    }
  }

  for (int i = 0; i < names.size(); ++i) {
    ss << "{\"n\": \"" << names[i] << "\", \"c\":" << comps[i] << "}";
    if (i < (names.size() - 1)) {
      ss << ",";
    }
  }

  ss << "], \"coarsest_geom\": {\"prob_lo\": "
    << "["
    << geom[0].ProbLo()[0] << ", "
    << geom[0].ProbLo()[1] << ", "
    << geom[0].ProbLo()[2]
    << "], "
    << "\"prob_hi\": "
    << "["
    << geom[0].ProbHi()[0] << ", "
    << geom[0].ProbHi()[1] << ", "
    << geom[0].ProbHi()[2]
    << "], "
    << "\"dx\": "
    << "["
    << geom[0].CellSize()[0] << ", "
    << geom[0].CellSize()[1] << ", "
    << geom[0].CellSize()[2]
    << "]}";

  ss << "}}";

  header = ss.str();
}

void MFIXReadWrite::SetMPMDRoots ()
{
  int mpmd_proc = MPMD::MyProc();
  int rank_offset = mpmd_proc - ParallelDescriptor::MyProc();
  if (rank_offset == 0) { // First program
    m_mpmd_this_root = 0;
    m_mpmd_other_root = ParallelDescriptor::NProcs();
  } else { // Second program
    m_mpmd_this_root = rank_offset;
    m_mpmd_other_root = 0;
  }
}

void MFIXReadWrite::SendMPMDIntFlags (bool end_flag) const
{
  if (MPMD::MyProc() == m_mpmd_this_root) {
    int tag = 0; // arbitrary tag

    Vector<int> int_flags;
    int end = end_flag ? 1 : 0;
    int_flags.push_back(end);

    MPI_Send(int_flags.data(), int_flags.size(), MPI_INT,
        m_mpmd_other_root, tag, MPI_COMM_WORLD);
  }
}

void MFIXReadWrite::SendMPMDReals (Real time) const
{
  if (MPMD::MyProc() == m_mpmd_this_root) {
    int tag = 10; // arbitrary tag

    Vector<Real> reals;
    reals.push_back(time);

    MPI_Send(reals.data(), reals.size(), MPI_DOUBLE, m_mpmd_other_root,
        tag, MPI_COMM_WORLD);
  }
}

void MFIXReadWrite::SendMPMDMultiFabs () const
{
  for (int lev = 0; lev < nlev; ++lev) {
    for (const auto& name: m_mpmd_mf_names) {
      if (name == "vel_g") {
        m_mpmd_copiers[lev]->send(*(leveldata_const().vel_const(lev)),
            0, AMREX_SPACEDIM);
      } else if (name == "ep_g") {
        m_mpmd_copiers[lev]->send( *(leveldata_const().epf_const(lev)), 0, 1);
      } else if (name == "T_g") {
        m_mpmd_copiers[lev]->send( *(leveldata_const().T_const(lev)), 0, 1);
      } else if (name == "X_gk") {
        for (int n = 0; n < fluid.nspecies(); ++n) {
          m_mpmd_copiers[lev]->send( *(leveldata_const().X_const(lev)), n, 1);
        }
      }
    }
  }
}

void MFIXReadWrite::MPMDInit ()
{
  Print() << "  Initializing MPMD\n";

  // TODO: Currently only works on 1 level
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nlev == 1, "MPMD currently only works with 1 level");

  // Variable check
  for (const auto& name: m_mpmd_mf_names) {
    if ((name == "T_g") && (!fluid.solve_enthalpy())) {
        Abort("Cannot send T_g via mpmd without solving enthalpy.");
    } else if ((name == "X_gk") && (!fluid.solve_species())) {
        Abort("Cannot send X_gk via mpmd without solving species.");
    }
  }

  // Define the copiers
  m_mpmd_copiers.resize(nlev);

  for (int lev = 0; lev < nlev; ++lev) {
    m_mpmd_copiers[lev] = std::make_unique<MPMD::Copier>(grids[lev],
          dmap[lev],true);
  }

  // Send the header info as a json string from this root
  // to the other root
  SetMPMDRoots();

  std::string header_json;
  ConstructMPMDHeader(header_json);

  if (MPMD::MyProc() == m_mpmd_this_root) {
    MPI_Send(header_json.c_str(), header_json.size(),
        MPI_CHAR, m_mpmd_other_root, 0, MPI_COMM_WORLD);
  }

  // Send the static multifabs once
  for (int lev = 0; lev < nlev; ++lev) {
    for (const auto& name: m_mpmd_static_mf_names) {
      if (name == "volfrac") {
        m_mpmd_copiers[lev]->send(ebfactory[lev]->getVolFrac(), 0, 1);
      } else if (name == "centroid") {
        m_mpmd_copiers[lev]->send(ebfactory[lev]->getCentroid().ToMultiFab(0., 0.),
            0, 3);
      }
    }
  }
}

void MFIXReadWrite::MPMDSend (Real time) const
{
  Print() << "  Sending MPMD data at time " << time << "\n";
  SendMPMDIntFlags(false);
  SendMPMDReals(time);
  SendMPMDMultiFabs();
}

void MFIXReadWrite::MPMDFinal () const
{
  Print() << "  Finalizing MPMD\n";
  SendMPMDIntFlags(true);
}
