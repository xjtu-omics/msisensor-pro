// System header files
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string.h>

#include "cmds.h"

#ifndef VERSION
#define VERSION "dev"
#endif

static void
print_version (void)
{
  const char *v = VERSION;
  std::cout << "msisensor-pro " << v << "\n";
}

static bool
is_version_flag (const char *arg)
{
  return strcmp (arg, "-v") == 0 || strcmp (arg, "-V") == 0
         || strcmp (arg, "--version") == 0 || strcmp (arg, "version") == 0;
}

int
usage (std::ostream &os = std::cerr)
{
  os << "\n\n"
     << " msisensor-pro: Microsatellite Instability (MSI) detection using "
        "high-throughput sequencing data. \n"
     << "         (Support tumor-normal paired samples and tumor-only "
        "samples) "
        "\n\n"
     << "Version: ";
  if (VERSION[0] == 'v' || VERSION[0] == 'V')
    {
      os << VERSION;
    }
  else
    {
      os << "v" << VERSION;
    }
  os << "\n\n"
     << "Usage: msisensor-pro <command> [options]\n\n"
     << "Commands:\n\n"
     << "\t scan\n"
     << "\t   scan the reference genome to get microsatellites information\n\n"
     << "\t msi\n"
     << "\t   evaluate MSI using paired tumor-normal sequencing data\n\n"
     << "\t baseline\n"
     << "\t   build baseline for tumor only detection\n\n"
     << "\t pro\n"
     << "\t   evaluate MSI using single (tumor) sample sequencing data \n\n"
     << "Options / Commands:\n\n"
     << "\t -v, -V, --version, version\n"
     << "\t   show version information\n\n"
     << "\t -h, --help\n"
     << "\t   show this help message\n\n"
     << "    If you have any questions, please open an issue on GitHub "
        "(https://github.com/xjtu-omics/msisensor-pro) or contact with Peng "
        "Jia (pengjia@xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn) .\n\n";

  return 1;
}

int
main (int argc, char **argv)
{
  try
    {
      if (argc < 2)
        {
          usage (std::cerr);
          return 1;
        }

      if (is_version_flag (argv[1]))
        {
          print_version ();
          return 0;
        }

      if (strcmp (argv[1], "-h") == 0 || strcmp (argv[1], "--help") == 0)
        {
          usage (std::cout);
          return 0;
        }

      if (argc >= 3 && is_version_flag (argv[2]))
        {
          print_version ();
          return 0;
        }

      if (strcmp (argv[1], "scan") == 0)
        {
          HomoAndMicrosateScan (argc - 1, argv + 1);
          return 0;
        }
      else if (strcmp (argv[1], "baseline") == 0)
        {
          TrainMsiP (argc - 1, argv + 1);
          return 0;
        }
      else if (strcmp (argv[1], "msi") == 0)
        {
          HomoAndMicrosateDisMsi (argc - 1, argv + 1);
          return 0;
        }
      else if (strcmp (argv[1], "pro") == 0)
        {
          HomoAndMicrosateDisMsiPro (argc - 1, argv + 1);
          return 0;
        }
      else
        {
          std::cerr << "ERROR: unrecognized command " << argv[1] << "\n";
          return 2;
        }
    }
  catch (const char *e)
    {
      std::cerr << e << std::endl;
      return 3;
    }

  return 0;
}
