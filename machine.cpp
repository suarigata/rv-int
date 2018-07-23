#include <elfio/elfio.hpp>
#include <machine.hpp>

#include <cstring>
#include <sstream>
#include <set>
#include <iomanip>

using namespace dbt;

//#define MAX_ARGUMENT_SIZE 1024 * 1024 /* 1mb */

union HalfUn {
	char asC_[2];
	uint16_t asH_;
};

void copystr(char* Target, const char* Source, uint32_t Size) {
  for (uint32_t i = 0; i < Size; ++i)
    Target[i] = Source[i];
}

void Machine::setCodeMemory(uint32_t StartAddress, uint32_t Size, const char* CodeBuffer) {
  CodeMemOffset = StartAddress;
  CodeMemory = uptr<Word[]>(new Word[Size]);
  CodeMemLimit = Size + CodeMemOffset;
  for (uint32_t i = 0; i < Size; i++) {
    Word Bytes = {CodeBuffer[i], CodeBuffer[i+1], CodeBuffer[i+2], CodeBuffer[i+3]};
    CodeMemory[i] = Bytes;
  }
}

void Machine::allocDataMemory(uint32_t Offset, uint32_t TotalSize) {
  DataMemTotalSize = TotalSize;
  DataMemOffset = Offset;
  DataMemLimit = Offset + TotalSize;
  DataMemory = std::unique_ptr<char[]>(new char[TotalSize]);
}

void Machine::addDataMemory(uint32_t StartAddress, uint32_t Size, const char* DataBuffer) {
  uint32_t Offset = StartAddress - DataMemOffset;
  DataMemLimit += Size;               //Ops, allocated memory stills the same, no more allocation is done and DataMemLimit is updated!
  copystr(DataMemory.get() + Offset, DataBuffer, Size);
}

int Machine::setCommandLineArguments(std::string parameters) {
  unsigned int sp = getRegister(2), totalSize=0, offset;

  std::istringstream iss(parameters);
  std::vector<std::string> argv(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
  argv.insert(argv.cbegin(), BinPath);

  for (auto argument : argv)
    totalSize += argument.length()+1;

  offset = DataMemTotalSize-totalSize-1;
  setMemValueAt(sp, (uint32_t) argv.size());

  for(auto argument : argv) {
    sp += 4;                                                          //Subtract stack pointer
    unsigned argSize = argument.length()+1;                           //Argument size
    copystr(DataMemory.get() + offset, argument.c_str(), argSize);    //Put argument in sp+4+size(arg[0..])=offset
    setMemValueAt(sp, (uint32_t) offset+DataMemOffset);               //Put offset in sp
    offset += argSize;                                                //Increment offset by argument Size
  }

  setMemValueAt(sp+4, 0);
  copystr(DataMemory.get() + offset, "\0", 1);
  return 0;
}

uint32_t Machine::getPC() {
  return PC;
}

uint32_t Machine::getLastPC() {
  return LastPC;
}

void Machine::incPC() {
  LastPC=PC;
  PC += 4;
}

void Machine::setPC(uint32_t NewPC) {
  LastPC = PC;
  PC = NewPC;
  #ifdef PRINTREG
  std::cerr << "LastPC= " << LastPC << "; NewPC= " << PC << ";  ";
  #endif
}

dbt::Word Machine::getInstAt(uint32_t Addr) {
  return CodeMemory[Addr - CodeMemOffset];
}

dbt::Word Machine::getInstAtPC() {
  return getInstAt(PC);
}

dbt::Word Machine::getNextInst() {
  ++PC;
  return getInstAtPC();
}

void Machine::setMemByteAt(uint32_t Addr, uint8_t Value) {
  uint32_t CorrectAddr = Addr - DataMemOffset;
  DataMemory[CorrectAddr] = Value;
}

uint8_t Machine::getMemByteAt(uint32_t Addr) {
  uint32_t CorrectAddr = Addr - DataMemOffset;
  return DataMemory[CorrectAddr];
}

uint16_t Machine::getMemHalfAt(uint32_t Addr) {
  uint32_t CorrectAddr = Addr - DataMemOffset;
  HalfUn Half = {DataMemory[CorrectAddr], DataMemory[CorrectAddr+1]};
  return Half.asH_;
}

dbt::Word Machine::getMemValueAt(uint32_t Addr) {
  uint32_t CorrectAddr = Addr - DataMemOffset;
  if(Addr<DataMemOffset){
    std::cout << "antes\n";
  }
  if(Addr>DataMemLimit){
    std::cout << "depois\n";
  }
  Word Bytes;
  Bytes.asI_ = *((uint32_t*)(DataMemory.get() + CorrectAddr));
  return Bytes;
}

dbt::DWord Machine::getDMemValueAt(uint32_t Addr) {
  uint64_t CorrectAddr = Addr - DataMemOffset;
  DWord Bytes;
  //Bytes.asI32_[0] = *((uint32_t*)(DataMemory.get() + CorrectAddr)); // TODO checar
  //Bytes.asI32_[1] = *((uint32_t*)(DataMemory.get() + CorrectAddr + 4));
  Bytes.asI_ = *((uint64_t*)(DataMemory.get() + CorrectAddr));
  return Bytes;
}

void Machine::setMemValueAt(uint32_t Addr, uint32_t Value) {
  uint32_t CorrectAddr = Addr - DataMemOffset;
  *((uint32_t*)(DataMemory.get() + CorrectAddr)) = Value;
}

void Machine::setDMemValueAt(uint32_t Addr, uint64_t Value) {
  uint64_t CorrectAddr = Addr - DataMemOffset;
  *((uint64_t*)(DataMemory.get() + CorrectAddr)) = Value;
}

uint32_t Machine::getNumInst() {
  return (CodeMemLimit - CodeMemOffset)/4;
}

uint32_t Machine::getCodeStartAddrs() {
  return CodeMemOffset;
}

uint32_t Machine::getCodeEndAddrs() {
  return CodeMemLimit;
}

uint32_t Machine::getDataMemOffset() {
  return DataMemOffset;
}

int32_t Machine::getRegister(uint16_t R) {
  #ifdef PRINTREG
  std::cerr << "GET R[" << std::dec << R << "]: " << std::hex << Register[R] << ";  ";
  #endif
  return R ? Register[R] : 0;
}

float Machine::getFloatRegister(uint16_t R) {
  return ((float*) Register)[R + 32]; // TODO TODO TODO deve ser 32 aqui!!!!
}

double Machine::getDoubleRegister(uint16_t R) {
  return ((double*)Register)[R + 65];
}

void Machine::setRegister(uint16_t R, int32_t V) {
    #ifdef PRINTREG
    std::cerr << "R[" << std::dec << R << "] = " << std::hex << V << ";  ";
    #endif
    Register[R] = V;
}

void Machine::setFloatRegister(uint16_t R, float V) {
  ((float*)Register)[R + 32] = V; // TODO era 66
}

void Machine::setDoubleRegister(uint16_t R, double V) {
  ((double*)Register)[R + 65] = V;
}

int32_t* Machine::getRegisterPtr() {
  return Register;
}

char* Machine::getByteMemoryPtr() {
  return DataMemory.get();
}

uint32_t* Machine::getMemoryPtr() {
  return (uint32_t*) DataMemory.get();
}

bool Machine::isOnNativeExecution() {
  return OnNativeExecution;
}

uint32_t Machine::getRegionBeingExecuted() {
  return RegionBeingExecuted;
}

void Machine::setOnNativeExecution(uint32_t EntryRegionAddrs) {
  OnNativeExecution   = true;
  RegionBeingExecuted = EntryRegionAddrs;
}

void Machine::setOffNativeExecution() {
  OnNativeExecution   = false;

}

uint32_t Machine::findMethod(uint32_t Addr) {
  for (auto Method : Symbolls) 
    if (Method.first < Addr && Method.second.second > Addr) 
      return Method.first;
  return 0;
}

bool Machine::isMethodEntry(uint32_t Addr) {
  return Symbolls.count(Addr) != 0;
}

uint32_t Machine::getMethodEnd(uint32_t Addr) {
  return Symbolls[Addr].second;
}

std::string Machine::getMethodName(uint32_t Addr) {
  return Symbolls[Addr].first;
}

std::vector<uint32_t> Machine::getVectorOfMethodEntries() {
  std::vector<uint32_t> R;
  for (auto KV : Symbolls)
    R.push_back(KV.first);
  return R;
}

using namespace ELFIO;

void Machine::reset() {
  loadELF(BinPath);
}

int Machine::loadELF(const std::string ElfPath) {
  BinPath = ElfPath;

  elfio reader;

  if (!reader.load(ElfPath))
    return 0;

  Elf_Half sec_num = reader.sections.size();

  uint32_t TotalDataSize = 0;
  uint32_t AddressOffset = 0;
  bool Started = false;
  bool First = false;
  for (int i = 0; i < sec_num; ++i) {
    section* psec = reader.sections[i];

    if (Started && (psec->get_flags() & 0x2) != 0) {
      TotalDataSize += psec->get_size();
      if (!First) {
        AddressOffset = psec->get_address();
        First = true;
      }
    }

    if (psec->get_name() == ".text")
      Started = true;
  }

  allocDataMemory(AddressOffset, (TotalDataSize + stackSize + heapSize) + (4 - (TotalDataSize + stackSize + heapSize) % 4));

  std::unordered_map<uint32_t, std::string> SymbolNames;
  std::set<uint32_t> SymbolStartAddresses;

  Started = false;
  for (int i = 0; i < sec_num; ++i) {
    section* psec = reader.sections[i];
    if (Started && (psec->get_flags() & 0x2) != 0 && psec->get_data() != nullptr) {
      addDataMemory(psec->get_address(), psec->get_size(), psec->get_data());
    }

    if (psec->get_name() == ".text") {
      setCodeMemory(psec->get_address(), psec->get_size(),  psec->get_data());
      SymbolStartAddresses.insert(psec->get_address() + psec->get_size());
      Started = true;
    }

    if (psec->get_name() == ".symtab") {
      const symbol_section_accessor symbols(reader, psec);
      std::string   name = "";
      Elf64_Addr    value = 0;
      Elf_Xword     size;
      unsigned char bind;
      unsigned char type = 0;
      Elf_Half      section_index;
      unsigned char other;
      for ( unsigned int j = 0; j < symbols.get_symbols_num(); ++j ) {
        symbols.get_symbol( j, name, value, size, bind, type, section_index, other );
        if (type == 0 && name != "" && value != 0 && value < CodeMemLimit) {
          SymbolStartAddresses.insert(value);
          SymbolNames[value] = name;
        }
      }
    }
  }

  for (auto I = SymbolStartAddresses.begin(); I != SymbolStartAddresses.end(); ++I)
    Symbolls[*I] = {SymbolNames[*I], *SymbolStartAddresses.upper_bound(*I)};

  for (int i = 0; i < 258; i++)
    Register[i] = 0;

  uint32_t StackAddr = DataMemLimit-stackSize/4;
  setRegister(2, StackAddr + (4 - StackAddr%4)); // StackPointer
  setRegister(30, StackAddr + (4 - StackAddr%4)); // StackPointer TODO que isso?

  setPC(reader.get_entry());

  return 1;
}

//#ifdef DEBUG
void Machine::dumpRegisters(void) {
  std::cerr << std::endl << "PC: " << std::hex << PC << "; \n";
  std::cerr << "(int32_t [258]) = {" << std::endl;

  for (int i = 0; i<258; ++i) {
    std::cerr << "  [" << std::dec << i << "] = " << "0x" << std::setw(8) << std::setfill('0') << std::hex <<  Register[i] << std::endl;
  }

  std::cerr << "}" << std::endl;
}
//#endif
