#include <arglib/arglib.hpp>
#include <interpreter.hpp>
#include <syscall.hpp>

#include <iostream>
#include <memory>

#include <signal.h>

clarg::argString RFTFlag("-rft", "Region Formation Technique (net)", "net");
clarg::argBool   InterpreterFlag("-interpret",  "Only interpret.");
clarg::argString BinaryFlag("-bin",  "path to the binary which will should be emulated.", "");
clarg::argBool   VerboseFlag("-v",  "display the compiled regions");
clarg::argBool   HelpFlag("-h",  "display the help message");
clarg::argInt    RegionLimitSize("-l", "region size limit", 400);
clarg::argString ArgumentsFlag("-args", "Pass Parameters to binary file (as string)", "");

void usage(char* PrgName) {
  cout << "Version: 0.0.1 (03-01-2017)\n\n";

  cout << "Usage: " << PrgName << 
    " [-rft net] [-interpreter] -bin PathToBinary\n\n";

  cout << "DESCRIPTION:\n";
  cout << "This program implements the OpenISA DBT (Dynamic Binary Translator)\n" <<
    "Institute of Computing, 2017.\n\n";

  cout << "ARGUMENTS:\n";
  clarg::arguments_descriptions(cout, "  ", "\n");
}

int validateArguments() {
  if (InterpreterFlag.was_set() && RFTFlag.was_set()) {
    cerr << "You can't use a RFT when you are only interpreting!\n";
    return 1;
  }

  if (!BinaryFlag.was_set()) {
    cerr << "You must set the path of the binary which will be emulated!\n";
    return 1;
  }

  return 0;
}

dbt::Machine M;

void  sigHandler(int sig) {
  std::cerr << "Segfault (" << sig << ") while emulating at PC: " << std::hex << M.getPC() << "\n";

  //TODO: Implement a M.dump();
  exit(1);
}

int main(int argc, char** argv) {
  signal(SIGSEGV, sigHandler);

  // Parse the arguments
  if (clarg::parse_arguments(argc, argv)) {
    cerr << "Error when parsing the arguments!" << endl;
    return 1;
  }

  if (HelpFlag.get_value() == true) {
    usage(argv[0]);
    return 1;
  }

  if (validateArguments())
    return 1;

  int loadStatus = M.loadELF(BinaryFlag.get_value());

  if (!loadStatus) {
    std::cout << "Can't find or process ELF file " << argv[1] << std::endl;
    return 2;
  }

  std::unique_ptr<dbt::SyscallManager> SyscallM;
  SyscallM = std::make_unique<dbt::LinuxSyscallManager>();

  if (M.setCommandLineArguments(ArgumentsFlag.get_value()) < 0)                                                                                  
    exit(1);  

  dbt::ITDInterpreter I(*SyscallM.get());
  std::cout << "Starting execution:\n";
  I.executeAll(M);

  return SyscallM->getExitStatus();
}

/*
 * TODO:
 *  - Add all instructions need by printf (commit)
 *  - Correct bugs in at least 6 tests (commit)
 *  - Add new benchmarks (commit)
 *  - Make improvements on the code (commit)
 *  - Make improvements on the performance (commit)
 *  ---------------------------------------------------- Until: 24 May
*/
