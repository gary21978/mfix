#include <fml_plotfile.H>

#include <AMReX_Print.H>

using namespace amrex;

fml_plotfile::
fml_plotfile ( std::string a_plot_file_name )
  : PlotFileData( a_plot_file_name )
  , m_pf_name(a_plot_file_name)
  , m_has_particles(false)
{

// Process particle data

  std::string hdr_file_name(m_pf_name + "/particles/Header");

  std::ifstream file(hdr_file_name.c_str());
  if ( file.is_open() ) {

     m_has_particles = true;

    // This is where we consume the data in the particle header file
    file >> *this;
    file.close();
  }

#ifdef BL_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}


std::istream& operator>> (std::istream& a_stream, fml_plotfile& a_header) {

  std::string pf_version;

  a_stream >> pf_version;
  a_stream >> a_header.m_pf_dim;

  a_stream >> a_header.m_pf_size_r_extra;
  for (int i = 0; i < a_header.m_pf_size_r_extra; ++i) {
    std::string comp;
    a_stream >> comp;
    a_header.m_pf_variables_r.push_back(comp);
  }

  a_stream >> a_header.m_pf_size_i;
  for (int i = 0; i < a_header.m_pf_size_i; ++i) {
    std::string comp;
    a_stream >> comp;
    a_header.m_pf_variables_i.push_back(comp);
  }

  a_stream >> a_header.m_is_checkpoint;
  a_stream >> a_header.m_nparticles;

  return a_stream;
}
