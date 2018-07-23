#include <interpreter.hpp>

#include <cmath>
#include <iostream>

using namespace dbt;
using namespace dbt::RVDecoder;

#ifdef COLOR
#define COLOR_CYAN "\033[1;35m"
#define COLOR_NONE "\033[0m"
#else
#define COLOR_CYAN ""
#define COLOR_NONE ""
#endif

//#define PRINTINST

uint64_t instacc = 0;
#ifdef PRINTINST
#include <RVPrinter.hpp>
//#define DEBUG_PRINT(Addr, Inst) std::cerr << RVPrinter::getString(Inst) << "\n";
//#define DEBUG_PRINT(Addr, Inst) std::cerr << /*std::dec << (++instacc) <<" -- "<<*/ std::hex << Addr << "\t" << RVPrinter::getString(Inst) << std::dec << "\n";
//#define DEBUG_PRINT(Addr, Inst) std::cerr << "\n" << COLOR_CYAN <<  std::dec << (++instacc) <<" -- "<< std::hex << Addr << "\t" << RVPrinter::getString(Inst) << std::dec << COLOR_NONE << "\t";
//#define DEBUG_PRINT(Addr, Inst) std::cerr << "\n" << COLOR_CYAN <<  std::dec << (++instacc) <<" -- "<< std::hex << Addr << "\t" << RVPrinter::getString(Inst) << std::dec << COLOR_NONE << "\t";
#else
#define DEBUG_PRINT(Addr, Inst)\
  {\
    std::cout << "Addr:" << std::hex << Addr << "\tImm: " << Inst.Imm\
              << std::hex << " RD " << Inst.RD\
              << std::hex << "\tRS1 " << Inst.RS1\
              << std::hex << "\tRS2 " << Inst.RS2\
              << std::hex << "\tRS3 " << Inst.RS3\
              << std::hex << "\tRM " << Inst.RM\
              << std::hex << " RD: " << M.getRegister(Inst.RD)  << "\n";\
  }
//TODO*/
#endif

#define SET_DISPACH(Addrs, Label, Offset)\
  case Label:\
    setDispatchValue(Addrs, static_cast<int*>(Offset));\
    break;

#define GOTO_NEXT\
    DEBUG_PRINT(M.getLastPC(), I)\
    goto *getDispatchValue(M.getPC());

#define IMPLEMENT(Label, Code)\
  Label:\
    I = getDecodedInst(M.getPC());\
    {\
       Code\
    }\
    M.incPC();\
    GOTO_NEXT

#define IMPLEMENT_JMP(Label, Code)\
  Label:\
    /*M.dumpRegisters();*/\
    I = getDecodedInst(M.getPC());\
    {\
      Code\
    }\
    GOTO_NEXT

#define IMPLEMENT_BR(Label, Code)\
  Label:\
    /*M.dumpRegisters();*/\
    I = getDecodedInst(M.getPC());\
    {\
      Code\
    }\
    M.incPC();\
    GOTO_NEXT

//#define rotate_right(x, n) (((x) >> (n)) | ((x) << ((sizeof(x) * 8) - (n))))
static inline uint32_t rotate_right(uint32_t input, uint32_t shiftamount) {
  return (((uint32_t)input) >> shiftamount) |
         (((uint32_t)input) << (32 - shiftamount));
}

static inline int32_t extend_signal(int32_t value, char type){
  switch(type){
    case 'I':
    case 'S': return (value & 0x800?0xFFFFF000|value:value);
    case 'B': return (value & 0x1000?0xFFFFE000|value:value);
    case 'J': return (value & 0x100000?0xFFE00000|value:value);
  }
  return 0;
}

typedef union DWordBit {
  double asF;
  uint64_t asI;
} DWordBit;

typedef union WordBit {
  float asF;
  int32_t asI;
} WordBit;

bool isnan(double x) { return x != x; }
bool isnan(float x)  { return x != x; }

bool ITDInterpreter::isAddrsContainedIn(uint32_t StartAddrs, uint32_t EndAddrs) {
  return !(StartAddrs < LastStartAddrs || EndAddrs > LastEndAddrs);
}

inline void* ITDInterpreter::getDispatchValue(uint32_t Addrs) {
  return (void*)DispatchValues[(Addrs-LastStartAddrs)/4];
}

void ITDInterpreter::setDispatchValue(uint32_t Addrs, int* Target) {
  DispatchValues[(Addrs-LastStartAddrs)/4] = Target;
}

inline RVInst ITDInterpreter::getDecodedInst(uint32_t Addrs) {
  return DecodedInsts[(Addrs-LastStartAddrs)/4];
}

void ITDInterpreter::setDecodedInst(uint32_t Addrs, RVInst DI) {
  DecodedInsts[(Addrs-LastStartAddrs)/4] = DI;
}

void ITDInterpreter::dispatch(Machine& M, uint32_t StartAddrs, uint32_t EndAddrs) {
  for (uint32_t Addrs = StartAddrs; Addrs < EndAddrs; Addrs+=4) {
    Word W = M.getInstAt(Addrs);
    RVInst I = decode(W.asI_);
    switch(I.Type) {
      SET_DISPACH(Addrs, LUI, &&lui);
      SET_DISPACH(Addrs, AUIPC, &&auipc);
      SET_DISPACH(Addrs, JAL, &&jal);
      SET_DISPACH(Addrs, JALR, &&jalr);
      SET_DISPACH(Addrs, BEQ, &&beq);
      SET_DISPACH(Addrs, BNE, &&bne);
      SET_DISPACH(Addrs, BLT, &&blt);
      SET_DISPACH(Addrs, BGE, &&bge);
      SET_DISPACH(Addrs, BLTU, &&bltu);
      SET_DISPACH(Addrs, BGEU, &&bgeu);
      SET_DISPACH(Addrs, LB, &&lb);
      SET_DISPACH(Addrs, LH, &&lh);
      SET_DISPACH(Addrs, LW, &&lw);
      SET_DISPACH(Addrs, LBU, &&lbu);
      SET_DISPACH(Addrs, LHU, &&lhu);
      SET_DISPACH(Addrs, SB, &&sb);
      SET_DISPACH(Addrs, SH, &&sh);
      SET_DISPACH(Addrs, SW, &&sw);
      SET_DISPACH(Addrs, ADDI, &&addi);
      SET_DISPACH(Addrs, SLTI, &&slti);
      SET_DISPACH(Addrs, SLTIU, &&sltiu);
      SET_DISPACH(Addrs, XORI, &&xori);
      SET_DISPACH(Addrs, ORI, &&ori);
      SET_DISPACH(Addrs, ANDI, &&andi);
      SET_DISPACH(Addrs, SLLI, &&slli);
      SET_DISPACH(Addrs, SRLI, &&srli);
      SET_DISPACH(Addrs, SRAI, &&srai);
      SET_DISPACH(Addrs, ADD, &&add);
      SET_DISPACH(Addrs, SUB, &&sub);
      SET_DISPACH(Addrs, SLL, &&sll);
      SET_DISPACH(Addrs, SLT, &&slt);
      SET_DISPACH(Addrs, SLTU, &&sltu);
      SET_DISPACH(Addrs, XOR, &&_xor);
      SET_DISPACH(Addrs, SRL, &&srl);
      SET_DISPACH(Addrs, SRA, &&sra);
      SET_DISPACH(Addrs, OR, &&_or);
      SET_DISPACH(Addrs, AND, &&_and);
      SET_DISPACH(Addrs, FENCE, &&fence);
      SET_DISPACH(Addrs, FENCE_I, &&fence_i);
      SET_DISPACH(Addrs, ECALL, &&ecall);
      SET_DISPACH(Addrs, EBREAK, &&ebreak);
      SET_DISPACH(Addrs, CSRRW, &&csrrw);
      SET_DISPACH(Addrs, CSRRS, &&csrrs);
      SET_DISPACH(Addrs, CSRRC, &&csrrc);
      SET_DISPACH(Addrs, CSRRWI, &&csrrwi);
      SET_DISPACH(Addrs, CSRRSI, &&csrrsi);
      SET_DISPACH(Addrs, CSRRCI, &&csrrci);
      SET_DISPACH(Addrs, MUL, &&mul);
      SET_DISPACH(Addrs, MULH, &&mulh);
      SET_DISPACH(Addrs, MULHSU, &&mulhsu);
      SET_DISPACH(Addrs, MULHU, &&mulhu);
      SET_DISPACH(Addrs, DIV, &&div);
      SET_DISPACH(Addrs, DIVU, &&divu);
      SET_DISPACH(Addrs, REM, &&rem);
      SET_DISPACH(Addrs, REMU, &&remu);
      SET_DISPACH(Addrs, FLW, &&flw);
      SET_DISPACH(Addrs, FSW, &&fsw);
      SET_DISPACH(Addrs, FMADD_S, &&fmadd_s);
      SET_DISPACH(Addrs, FMSUB_S, &&fmsub_s);
      SET_DISPACH(Addrs, FNMSUB_S, &&fnmsub_s);
      SET_DISPACH(Addrs, FNMADD_S, &&fnmadd_s);
      SET_DISPACH(Addrs, FADD_S, &&fadd_s);
      SET_DISPACH(Addrs, FSUB_S, &&fsub_s);
      SET_DISPACH(Addrs, FMUL_S, &&fmul_s);
      SET_DISPACH(Addrs, FDIV_S, &&fdiv_s);
      SET_DISPACH(Addrs, FSQRT_S, &&fsqrt_s);
      SET_DISPACH(Addrs, FSGNJ_S, &&fsgnj_s);
      SET_DISPACH(Addrs, FSGNJN_S, &&fsgnjn_s);
      SET_DISPACH(Addrs, FSGNJX_S, &&fsgnjx_s);
      SET_DISPACH(Addrs, FMIN_S, &&fmin_s);
      SET_DISPACH(Addrs, FMAX_S, &&fmax_s);
      SET_DISPACH(Addrs, FCVT_W_S, &&fcvt_w_s);
      SET_DISPACH(Addrs, FCVT_WU_S, &&fcvt_wu_s);
      SET_DISPACH(Addrs, FMV_X_W, &&fmv_x_w);
      SET_DISPACH(Addrs, FEQ_S, &&feq_s);
      SET_DISPACH(Addrs, FLT_S, &&flt_s);
      SET_DISPACH(Addrs, FLE_S, &&fle_s);
      SET_DISPACH(Addrs, FCLASS_S, &&fclass_s);
      SET_DISPACH(Addrs, FCVT_S_W, &&fcvt_s_w);
      SET_DISPACH(Addrs, FCVT_S_WU, &&fcvt_s_wu);
      SET_DISPACH(Addrs, FMV_W_X, &&fmv_w_x);
      SET_DISPACH(Addrs, FLD, &&fld);
      SET_DISPACH(Addrs, FSD, &&fsd);
      SET_DISPACH(Addrs, FMADD_D, &&fmadd_d);
      SET_DISPACH(Addrs, FMSUB_D, &&fmsub_d);
      SET_DISPACH(Addrs, FNMSUB_D, &&fnmsub_d);
      SET_DISPACH(Addrs, FNMADD_D, &&fnmadd_d);
      SET_DISPACH(Addrs, FADD_D, &&fadd_d);
      SET_DISPACH(Addrs, FSUB_D, &&fsub_d);
      SET_DISPACH(Addrs, FMUL_D, &&fmul_d);
      SET_DISPACH(Addrs, FDIV_D, &&fdiv_d);
      SET_DISPACH(Addrs, FSQRT_D, &&fsqrt_d);
      SET_DISPACH(Addrs, FSGNJ_D, &&fsgnj_d);
      SET_DISPACH(Addrs, FSGNJN_D, &&fsgnjn_d);
      SET_DISPACH(Addrs, FSGNJX_D, &&fsgnjx_d);
      SET_DISPACH(Addrs, FMIN_D, &&fmin_d);
      SET_DISPACH(Addrs, FMAX_D, &&fmax_d);
      SET_DISPACH(Addrs, FCVT_S_D, &&fcvt_s_d);
      SET_DISPACH(Addrs, FCVT_D_S, &&fcvt_d_s);
      SET_DISPACH(Addrs, FEQ_D, &&feq_d);
      SET_DISPACH(Addrs, FLT_D, &&flt_d);
      SET_DISPACH(Addrs, FLE_D, &&fle_d);
      SET_DISPACH(Addrs, FCLASS_D, &&fclass_d);
      SET_DISPACH(Addrs, FCVT_W_D, &&fcvt_w_d);
      SET_DISPACH(Addrs, FCVT_WU_D, &&fcvt_wu_d);
      SET_DISPACH(Addrs, FCVT_D_W, &&fcvt_d_w);
      SET_DISPACH(Addrs, FCVT_D_WU, &&fcvt_d_wu);
      case Null:
        exit(1);
    }
    setDecodedInst(Addrs, I);
  }

  // ---------------------------------------- Trampoline Zone ------------------------------------------ //

  RVInst I;
  GOTO_NEXT;

//  IMPLEMENT(nop, ); TODO sem nop?

  /**********************   Int Inst   **************************/

  IMPLEMENT(lui,
    M.setRegister(I.RD, I.Imm);
  );

  IMPLEMENT(auipc,
    M.setRegister(I.RD, M.getPC() + I.Imm);
  );

  IMPLEMENT_JMP(jal,
    M.setRegister(I.RD, M.getPC() + 4);
    M.setPC(M.getPC() + extend_signal(I.Imm, 'J'));
  );

  IMPLEMENT_JMP(jalr,
    uint32_t _PC=M.getPC();
    M.setPC((M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')) & 0xFFFFFFFE);
    M.setRegister(I.RD, _PC + 4);
  );

  IMPLEMENT_BR(beq,
    if(M.getRegister(I.RS1) == M.getRegister(I.RS2)){
      M.setPC(M.getPC() + extend_signal(I.Imm, 'B'));
      GOTO_NEXT;
    }
  );

  IMPLEMENT_BR(bne,
    if(M.getRegister(I.RS1) != M.getRegister(I.RS2)){
      M.setPC(M.getPC() + extend_signal(I.Imm, 'B'));
      GOTO_NEXT;
    }
  );

  IMPLEMENT_BR(blt,
    if(M.getRegister(I.RS1) < M.getRegister(I.RS2)){
      M.setPC(M.getPC() + extend_signal(I.Imm, 'B'));
      GOTO_NEXT;
    }
  );

  IMPLEMENT_BR(bge,
    if(M.getRegister(I.RS1) >= M.getRegister(I.RS2)){
      M.setPC(M.getPC() + extend_signal(I.Imm, 'B'));
      GOTO_NEXT;
    }
  );

  IMPLEMENT_BR(bltu,
    if((uint32_t)M.getRegister(I.RS1) < (uint32_t)M.getRegister(I.RS2)){
      M.setPC(M.getPC() + extend_signal(I.Imm, 'B'));
      GOTO_NEXT;
    }
  );

  IMPLEMENT_BR(bgeu,
    if((uint32_t)M.getRegister(I.RS1) >= (uint32_t)M.getRegister(I.RS2)){
      M.setPC(M.getPC() + extend_signal(I.Imm, 'B'));
      GOTO_NEXT;
    }
  );

  IMPLEMENT(lb,
    M.setRegister(I.RD, (int32_t)(int8_t)M.getMemByteAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')));
  );

  IMPLEMENT(lh,
    M.setRegister(I.RD, (int32_t)(int16_t)M.getMemHalfAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')));
  );

  IMPLEMENT(lw,
    M.setRegister(I.RD, M.getMemValueAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')).asI_);
  );

  IMPLEMENT(lbu,
    M.setRegister(I.RD, (uint32_t)(uint8_t)M.getMemByteAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')));
  );

  IMPLEMENT(lhu,
    M.setRegister(I.RD, (uint32_t)(uint16_t)M.getMemHalfAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')));
  );

  IMPLEMENT(sb,
    std::cout << "sbmem: " << M.getRegister(I.RS1) + extend_signal(I.Imm, 'S') << "\n";
    M.setMemByteAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'S'), (unsigned char) M.getRegister(I.RS2) & 0xFF);
  );

  IMPLEMENT(sh,
    uint16_t half = M.getRegister(I.RS2) & 0xFFFF;
    M.setMemByteAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'S')    , (half)      & 0xFF);
    M.setMemByteAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'S') + 1, (half >> 8) & 0xFF);
  );

  IMPLEMENT(sw,
    M.setMemValueAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'S'), M.getRegister(I.RS2));
  );

  IMPLEMENT(addi,
    M.setRegister(I.RD, M.getRegister(I.RS1) + extend_signal(I.Imm, 'I'));
  );

  IMPLEMENT(slti,
    M.setRegister(I.RD, (M.getRegister(I.RS1) < extend_signal(I.Imm, 'I') ? 1 : 0));
  );

  IMPLEMENT(sltiu,
    M.setRegister(I.RD, ((uint32_t)(M.getRegister(I.RS1)) < (uint32_t)(extend_signal(I.Imm, 'I')) ? 1 : 0));
  );

  IMPLEMENT(xori,
    M.setRegister(I.RD, M.getRegister(I.RS1) ^ extend_signal(I.Imm, 'I'));
  );

  IMPLEMENT(ori,
    M.setRegister(I.RD, M.getRegister(I.RS1) | extend_signal(I.Imm, 'I'));
  );

  IMPLEMENT(andi,
    M.setRegister(I.RD, M.getRegister(I.RS1) & extend_signal(I.Imm, 'I'));
  );

  IMPLEMENT(slli,
    M.setRegister(I.RD, M.getRegister(I.RS1) << (I.Imm & 0x3F));
  );

  IMPLEMENT(srli,
    M.setRegister(I.RD, M.getRegister(I.RS1) >> (I.Imm & 0x3F));
  );

  IMPLEMENT(srai,
    int32_t _REG=M.getRegister(I.RS1);
    M.setRegister(I.RD, ((_REG & 0x80000000) ? (_REG >> (I.Imm & 0x3F)) | (0xFFFFFFFF << (32-(I.Imm & 0x3F))) : _REG >> (I.Imm & 0x3F)));
  );

  IMPLEMENT(add,
    M.setRegister(I.RD, M.getRegister(I.RS1) + M.getRegister(I.RS2));
  );

  IMPLEMENT(sub,
    M.setRegister(I.RD, M.getRegister(I.RS1) - M.getRegister(I.RS2));
  );

  IMPLEMENT(sll,
    M.setRegister(I.RD, M.getRegister(I.RS1) << (M.getRegister(I.RS2) & 0x1F));
  );

  IMPLEMENT(slt,
    M.setRegister(I.RD, (M.getRegister(I.RS1) < M.getRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(sltu,
    M.setRegister(I.RD, ((uint32_t)(M.getRegister(I.RS1)) < (uint32_t)(M.getRegister(I.RS2)) ? 1 : 0));
  );

  IMPLEMENT(_xor,
    M.setRegister(I.RD, M.getRegister(I.RS1) ^ M.getRegister(I.RS2));
  );

  IMPLEMENT(srl,
      M.setRegister(I.RD, M.getRegister(I.RS1) >> (M.getRegister(I.RS2) & 0x1F));
  );

  IMPLEMENT(sra,
    int32_t _REG1=M.getRegister(I.RS1);
    int32_t _REG2=M.getRegister(I.RS2) & 0x1F;
    M.setRegister(I.RD, ((_REG1 & 0x80000000) ? (_REG1 >> _REG2) | (0xFFFFFFFF << (32-_REG2)) : _REG1 >> _REG2));
  );

  IMPLEMENT(_or,
    M.setRegister(I.RD, M.getRegister(I.RS1) | M.getRegister(I.RS2));
  );

  IMPLEMENT(_and,
    M.setRegister(I.RD, M.getRegister(I.RS1) & M.getRegister(I.RS2));
  );

  IMPLEMENT(fence,
    
  );

  IMPLEMENT(fence_i,
    
  );

  IMPLEMENT(ecall,
    if (SyscallM.processSyscall(M))
      return;
  );

  IMPLEMENT(ebreak,
    
  );

  IMPLEMENT(csrrw,
    
  );

  IMPLEMENT(csrrs,
    
  );

  IMPLEMENT(csrrc,
    
  );

  IMPLEMENT(csrrwi,
    
  );

  IMPLEMENT(csrrsi,
    
  );

  IMPLEMENT(csrrci,
    
  );

  IMPLEMENT(mul,
      M.setRegister(I.RD, M.getRegister(I.RS1) * M.getRegister(I.RS2));
  );

  IMPLEMENT(mulh,
      M.setRegister(I.RD, (int32_t)(((int64_t)(M.getRegister(I.RS1)) * (int64_t)(M.getRegister(I.RS2))) >> 32));
  );

  IMPLEMENT(mulhsu,
      M.setRegister(I.RD, (int32_t)(((int64_t)(int32_t)(M.getRegister(I.RS1)) * (uint64_t)(uint32_t)(M.getRegister(I.RS2))) >> 32));
  );

  IMPLEMENT(mulhu,
      M.setRegister(I.RD, (int32_t)(((uint64_t)(uint32_t)(M.getRegister(I.RS1)) * (uint64_t)(uint32_t)(M.getRegister(I.RS2))) >> 32));
  );

  IMPLEMENT(div,
      M.setRegister(I.RD, M.getRegister(I.RS1) / M.getRegister(I.RS2));
  );

  IMPLEMENT(divu,
      M.setRegister(I.RD, (int32_t)((uint32_t)(M.getRegister(I.RS1)) / (uint32_t)(M.getRegister(I.RS2))));
  );

  IMPLEMENT(rem,
      M.setRegister(I.RD, M.getRegister(I.RS1) % M.getRegister(I.RS2));
  );

  IMPLEMENT(remu,
      M.setRegister(I.RD, (int32_t)((uint32_t)(M.getRegister(I.RS1)) % (uint32_t)(M.getRegister(I.RS2))));
  );

  IMPLEMENT(flw,
    M.setFloatRegister(I.RD, M.getMemValueAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')).asF_);
  );

  IMPLEMENT(fsw,
    M.setMemValueAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'S'), M.getFloatRegister(I.RS2));
  );

  IMPLEMENT(fmadd_s,
    M.setFloatRegister(I.RD, M.getFloatRegister(I.RS1) * M.getFloatRegister(I.RS2) + M.getFloatRegister(I.RS3));
  );

  IMPLEMENT(fmsub_s,
    M.setFloatRegister(I.RD, M.getFloatRegister(I.RS1) * M.getFloatRegister(I.RS2) - M.getFloatRegister(I.RS3));
  );

  IMPLEMENT(fnmsub_s,
    M.setFloatRegister(I.RD, -M.getFloatRegister(I.RS1) * M.getFloatRegister(I.RS2) + M.getFloatRegister(I.RS3));
  );

  IMPLEMENT(fnmadd_s,
    M.setFloatRegister(I.RD, -M.getFloatRegister(I.RS1) * M.getFloatRegister(I.RS2) - M.getFloatRegister(I.RS3));
  );

  IMPLEMENT(fadd_s,
    M.setFloatRegister(I.RD, M.getFloatRegister(I.RS1) + M.getFloatRegister(I.RS2));
  );

  IMPLEMENT(fsub_s,
    M.setFloatRegister(I.RD, M.getFloatRegister(I.RS1) - M.getFloatRegister(I.RS2));
  );

  IMPLEMENT(fmul_s,
    M.setFloatRegister(I.RD, M.getFloatRegister(I.RS1) * M.getFloatRegister(I.RS2));
  );

  IMPLEMENT(fdiv_s,
    M.setFloatRegister(I.RD, M.getFloatRegister(I.RS1) / M.getFloatRegister(I.RS2));
  );

  IMPLEMENT(fsqrt_s,
    M.setFloatRegister(I.RD, sqrt(M.getFloatRegister(I.RS1)));
  );

  IMPLEMENT(fsgnj_s,
    float _RS1=M.getFloatRegister(I.RS1);
    float _RS2=M.getFloatRegister(I.RS2);
    M.setFloatRegister(I.RD, ((_RS1 < 0) == (_RS2 < 0)) ? _RS1 : -_RS1);
  );

  IMPLEMENT(fsgnjn_s,
    float _RS1=M.getFloatRegister(I.RS1);
    float _RS2=M.getFloatRegister(I.RS2);
    M.setFloatRegister(I.RD, ((_RS1 < 0) != (_RS2 < 0)) ? _RS1 : -_RS1);
  );

  IMPLEMENT(fsgnjx_s,
    float _RS1=M.getFloatRegister(I.RS1);
    float _RS2=M.getFloatRegister(I.RS2);
    M.setFloatRegister(I.RD, ((_RS1 < 0) == (_RS2 < 0)) ? (_RS1<0?-_RS1:_RS1) : (_RS1>0?-_RS1:_RS1));
  );

  IMPLEMENT(fmin_s,
    float _RS1=M.getFloatRegister(I.RS1);
    float _RS2=M.getFloatRegister(I.RS2);
    M.setFloatRegister(I.RD, _RS1 < _RS2 ? _RS1 : _RS2);
  );

  IMPLEMENT(fmax_s,
    float _RS1=M.getFloatRegister(I.RS1);
    float _RS2=M.getFloatRegister(I.RS2);
    M.setFloatRegister(I.RD, _RS1 > _RS2 ? _RS1 : _RS2);
  );

  IMPLEMENT(fcvt_w_s,
    M.setRegister(I.RD, (int32_t) M.getFloatRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fcvt_wu_s,
    M.setRegister(I.RD, (uint32_t) M.getFloatRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fmv_x_w,
    float _RS1=M.getFloatRegister(I.RS1);
    unsigned char* _CHR = reinterpret_cast<unsigned char*>(&_RS1);
    M.setRegister(I.RD, *reinterpret_cast<uint32_t*>(_CHR));
  );

  IMPLEMENT(feq_s,
    M.setRegister(I.RD, (M.getFloatRegister(I.RS1) == M.getFloatRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(flt_s,
    M.setRegister(I.RD, (M.getFloatRegister(I.RS1) < M.getFloatRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(fle_s,
    M.setRegister(I.RD, (M.getFloatRegister(I.RS1) <= M.getFloatRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(fclass_s,
    M.setRegister(I.RD, 0); // TODO
  );

  IMPLEMENT(fcvt_s_w,
    M.setFloatRegister(I.RD, (float) M.getRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fcvt_s_wu,
    M.setFloatRegister(I.RD, (float)(uint32_t) M.getRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fmv_w_x,
    //uint32_t _RS1=M.getRegister(I.RS1);
    Word RCD;
    RCD.asI_=M.getRegister(I.RS1);
    //unsigned char* _CHR = reinterpret_cast<unsigned char*>(&_RS1);
    //M.setFloatRegister(I.RD, *reinterpret_cast<float*>(_CHR));
    M.setFloatRegister(I.RD, RCD.asF_);
  );

  IMPLEMENT(fld,
    M.setDoubleRegister(I.RD, M.getDMemValueAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'I')).asD_);
  );

  IMPLEMENT(fsd,
    M.setDMemValueAt(M.getRegister(I.RS1) + extend_signal(I.Imm, 'S'), M.getDoubleRegister(I.RS2));
  );

  IMPLEMENT(fmadd_d,
    M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS1) * M.getDoubleRegister(I.RS2) + M.getDoubleRegister(I.RS3));
  );

  IMPLEMENT(fmsub_d,
    M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS1) * M.getDoubleRegister(I.RS2) - M.getDoubleRegister(I.RS3));
  );

  IMPLEMENT(fnmsub_d,
    M.setDoubleRegister(I.RD, -M.getDoubleRegister(I.RS1) * M.getDoubleRegister(I.RS2) + M.getDoubleRegister(I.RS3));
  );

  IMPLEMENT(fnmadd_d,
    M.setDoubleRegister(I.RD, -M.getDoubleRegister(I.RS1) * M.getDoubleRegister(I.RS2) - M.getDoubleRegister(I.RS3));
  );

  IMPLEMENT(fadd_d,
    M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS1) + M.getDoubleRegister(I.RS2));
  );

  IMPLEMENT(fsub_d,
    M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS1) - M.getDoubleRegister(I.RS2));
  );

  IMPLEMENT(fmul_d,
    M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS1) * M.getDoubleRegister(I.RS2));
  );

  IMPLEMENT(fdiv_d,
    M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS1) / M.getDoubleRegister(I.RS2));
  );

  IMPLEMENT(fsqrt_d,
    M.setDoubleRegister(I.RD, sqrt(M.getDoubleRegister(I.RS1)));
  );

  IMPLEMENT(fsgnj_d,
    double _RS1=M.getDoubleRegister(I.RS1);
    double _RS2=M.getDoubleRegister(I.RS2);
    M.setDoubleRegister(I.RD, ((_RS1 < 0) == (_RS2 < 0)) ? _RS1 : -_RS1);
  );

  IMPLEMENT(fsgnjn_d,
    double _RS1=M.getDoubleRegister(I.RS1);
    double _RS2=M.getDoubleRegister(I.RS2);
    M.setDoubleRegister(I.RD, ((_RS1 < 0) != (_RS2 < 0)) ? _RS1 : -_RS1);
  );

  IMPLEMENT(fsgnjx_d,
    double _RS1=M.getDoubleRegister(I.RS1);
    double _RS2=M.getDoubleRegister(I.RS2);
    M.setDoubleRegister(I.RD, ((_RS1 < 0) == (_RS2 < 0)) ? (_RS1<0?-_RS1:_RS1) : (_RS1>0?-_RS1:_RS1));
  );

  IMPLEMENT(fmin_d,
    double _RS1=M.getDoubleRegister(I.RS1);
    double _RS2=M.getDoubleRegister(I.RS2);
    M.setDoubleRegister(I.RD, _RS1 < _RS2 ? _RS1 : _RS2);
  );

  IMPLEMENT(fmax_d,
    double _RS1=M.getDoubleRegister(I.RS1);
    double _RS2=M.getDoubleRegister(I.RS2);
    M.setDoubleRegister(I.RD, _RS1 > _RS2 ? _RS1 : _RS2);
  );

  IMPLEMENT(fcvt_s_d,
    M.setFloatRegister(I.RD, (float) M.getDoubleRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fcvt_d_s,
    M.setDoubleRegister(I.RD, (double) M.getFloatRegister(I.RS1));
  );

  IMPLEMENT(feq_d,
    M.setRegister(I.RD, (M.getDoubleRegister(I.RS1) == M.getDoubleRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(flt_d,
    M.setRegister(I.RD, (M.getDoubleRegister(I.RS1) < M.getDoubleRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(fle_d,
    M.setRegister(I.RD, (M.getDoubleRegister(I.RS1) <= M.getDoubleRegister(I.RS2) ? 1 : 0));
  );

  IMPLEMENT(fclass_d,
    M.setRegister(I.RD, 0); // TODO
  );

  IMPLEMENT(fcvt_w_d,
    M.setRegister(I.RD, (int32_t) M.getDoubleRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fcvt_wu_d,
    M.setRegister(I.RD, (uint32_t) M.getDoubleRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fcvt_d_w,
    M.setDoubleRegister(I.RD, (double) M.getRegister(I.RS1)); // TODO
  );

  IMPLEMENT(fcvt_d_wu,
    M.setDoubleRegister(I.RD, (double)(uint32_t) M.getRegister(I.RS1)); // TODO
  );

  // --------------------------------------------------------------------------------------------------- //
}

void ITDInterpreter::execute(Machine& M, uint32_t StartAddrs, uint32_t EndAddrs) {
  if (DispatchValues.size() == 0 || !isAddrsContainedIn(StartAddrs, EndAddrs)) {
    DispatchValues.reserve((EndAddrs - StartAddrs)/4);
    DecodedInsts  .reserve((EndAddrs - StartAddrs)/4);
  }

  LastStartAddrs = StartAddrs;
  LastEndAddrs   = EndAddrs;

  dispatch(M, StartAddrs, EndAddrs);
}


