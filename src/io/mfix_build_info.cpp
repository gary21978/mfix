#include "AMReX_buildInfo.H"
#include "build_info.H"
#include <mfix.H>
#ifdef AMREX_USE_HYPRE
#include <HYPRE_config.h>
#endif


void writeBuildInfo ()
{
  std::string PrettyLine = std::string(78, '=') + "\n";
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // build information
  std::cout << PrettyLine;
  std::cout << " MFIX-Exa Build Information\n";
  std::cout << PrettyLine;

  std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
  std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
  std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
  std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  std::cout << "\n";

  std::cout << "COMP:          " << buildInfoGetComp() << "\n";
  std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  std::cout << "\n";

  std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  std::cout << "\n";

  std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";

  std::cout << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
  }

  std::cout << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  if (strlen(githash1) > 0) {
    std::cout << "MFIX-Exa git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    std::cout << "AMReX    git describe: " << githash2 << "\n";
  }

  std::cout << "AMReX-Hydro  git hash: " << HydroGitHash() << "\n";
#ifdef CSG_EB
  std::cout << "CSG-EB       git hash: " << CsgEbGitHash() << "\n";
#endif

  const char* buildgithash = buildInfoGetBuildGitHash();
  const char* buildgitname = buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0){
    std::cout << buildgitname << " git describe: " << buildgithash << "\n";
  }

  std::cout << "\n\n";
}

const char* HypreVersion () {
#ifdef HYPRE_DEVELOP_STRING
  return HYPRE_DEVELOP_STRING;
#elif defined HYPRE_RELEASE_VERSION
  return HYPRE_RELEASE_VERSION;
#else
  return "Unknown";
#endif
}
