#include <RVDecoder.hpp>
#include <machine.hpp>
#include <syscall.hpp>


namespace dbt {
  class Interpreter {

  protected:
    SyscallManager& SyscallM;

  public:
    Interpreter(SyscallManager& SM) : SyscallM(SM) {}

    virtual void execute(Machine&, uint32_t, uint32_t) = 0;

    void executeAll(Machine& M) {
      execute(M, M.getCodeStartAddrs(), M.getCodeEndAddrs());
    }
  };

  class ITDInterpreter : public Interpreter {
  private:
    uint32_t LastStartAddrs, LastEndAddrs;
    std::vector<int*> DispatchValues;
    std::vector<RVDecoder::OIInst> DecodedInsts;

    bool isAddrsContainedIn(uint32_t, uint32_t);

    void dispatch(Machine&, uint32_t, uint32_t);

    void* getDispatchValue(uint32_t);
    void setDispatchValue(uint32_t, int*);

    RVDecoder::OIInst getDecodedInst(uint32_t);
    void setDecodedInst(uint32_t, RVDecoder::OIInst);
  public:
    ITDInterpreter(SyscallManager& SM) : Interpreter(SM) {}

    void execute(Machine&, uint32_t, uint32_t);
  };
}
