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
#define DEBUG_PRINT(Addr, Inst)
#endif

#define SET_DISPACH(Addrs, Label, Offset)\
  case Label:\
    setDispatchValue(Addrs, static_cast<int*>(Offset));\
    break;

#define GOTO_NEXT\
    DEBUG_PRINT(M.getPC(), getDecodedInst(M.getPC()))\
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
    
  );

  IMPLEMENT(auipc,
    
  );

  IMPLEMENT(jal,
    
  );

  IMPLEMENT(jalr,
    
  );

  IMPLEMENT(beq,
    
  );

  IMPLEMENT(bne,
    
  );

  IMPLEMENT(blt,
    
  );

  IMPLEMENT(bge,
    
  );

  IMPLEMENT(bltu,
    
  );

  IMPLEMENT(bgeu,
    
  );

  IMPLEMENT(lb,
    
  );

  IMPLEMENT(lh,
    
  );

  IMPLEMENT(lw,
    
  );

  IMPLEMENT(lbu,
    
  );

  IMPLEMENT(lhu,
    
  );

  IMPLEMENT(sb,
    
  );

  IMPLEMENT(sh,
    
  );

  IMPLEMENT(sw,
    
  );

  IMPLEMENT(addi,
    
  );

  IMPLEMENT(slti,
    
  );

  IMPLEMENT(sltiu,
    
  );

  IMPLEMENT(xori,
    
  );

  IMPLEMENT(ori,
    
  );

  IMPLEMENT(andi,
    
  );

  IMPLEMENT(slli,
    
  );

  IMPLEMENT(srli,
    
  );

  IMPLEMENT(srai,
    
  );

  IMPLEMENT(add,
    
  );

  IMPLEMENT(sub,
    
  );

  IMPLEMENT(sll,
    
  );

  IMPLEMENT(slt,
    
  );

  IMPLEMENT(sltu,
    
  );

  IMPLEMENT(_xor,
    
  );

  IMPLEMENT(srl,
    
  );

  IMPLEMENT(sra,
    
  );

  IMPLEMENT(_or,
    
  );

  IMPLEMENT(_and,
    
  );

  IMPLEMENT(fence,
    
  );

  IMPLEMENT(fence_i,
    
  );

  IMPLEMENT(ecall,
    
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
    
  );

  IMPLEMENT(mulh,
    
  );

  IMPLEMENT(mulhsu,
    
  );

  IMPLEMENT(mulhu,
    
  );

  IMPLEMENT(div,
    
  );

  IMPLEMENT(divu,
    
  );

  IMPLEMENT(rem,
    
  );

  IMPLEMENT(remu,
    
  );

  IMPLEMENT(flw,
    
  );

  IMPLEMENT(fsw,
    
  );

  IMPLEMENT(fmadd_s,
    
  );

  IMPLEMENT(fmsub_s,
    
  );

  IMPLEMENT(fnmsub_s,
    
  );

  IMPLEMENT(fnmadd_s,
    
  );

  IMPLEMENT(fadd_s,
    
  );

  IMPLEMENT(fsub_s,
    
  );

  IMPLEMENT(fmul_s,
    
  );

  IMPLEMENT(fdiv_s,
    
  );

  IMPLEMENT(fsqrt_s,
    
  );

  IMPLEMENT(fsgnj_s,
    
  );

  IMPLEMENT(fsgnjn_s,
    
  );

  IMPLEMENT(fsgnjx_s,
    
  );

  IMPLEMENT(fmin_s,
    
  );

  IMPLEMENT(fmax_s,
    
  );

  IMPLEMENT(fcvt_w_s,
    
  );

  IMPLEMENT(fcvt_wu_s,
    
  );

  IMPLEMENT(fmv_x_w,
    
  );

  IMPLEMENT(feq_s,
    
  );

  IMPLEMENT(flt_s,
    
  );

  IMPLEMENT(fle_s,
    
  );

  IMPLEMENT(fclass_s,
    
  );

  IMPLEMENT(fcvt_s_w,
    
  );

  IMPLEMENT(fcvt_s_wu,
    
  );

  IMPLEMENT(fmv_w_x,
    
  );

  IMPLEMENT(fld,
    
  );

  IMPLEMENT(fsd,
    
  );

  IMPLEMENT(fmadd_d,
    
  );

  IMPLEMENT(fmsub_d,
    
  );

  IMPLEMENT(fnmsub_d,
    
  );

  IMPLEMENT(fnmadd_d,
    
  );

  IMPLEMENT(fadd_d,
    
  );

  IMPLEMENT(fsub_d,
    
  );

  IMPLEMENT(fmul_d,
    
  );

  IMPLEMENT(fdiv_d,
    
  );

  IMPLEMENT(fsqrt_d,
    
  );

  IMPLEMENT(fsgnj_d,
    
  );

  IMPLEMENT(fsgnjn_d,
    
  );

  IMPLEMENT(fsgnjx_d,
    
  );

  IMPLEMENT(fmin_d,
    
  );

  IMPLEMENT(fmax_d,
    
  );

  IMPLEMENT(fcvt_s_d,
    
  );

  IMPLEMENT(fcvt_d_s,
    
  );

  IMPLEMENT(feq_d,
    
  );

  IMPLEMENT(flt_d,
    
  );

  IMPLEMENT(fle_d,
    
  );

  IMPLEMENT(fclass_d,
    
  );

  IMPLEMENT(fcvt_w_d,
    
  );

  IMPLEMENT(fcvt_wu_d,
    
  );

  IMPLEMENT(fcvt_d_w,
    
  );

  IMPLEMENT(fcvt_d_wu,
    
  );

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO


  IMPLEMENT(add,

      M.setRegister(I.RD, M.getRegister(I.RS) + M.getRegister(I.RT));
    );

  IMPLEMENT(sub,
      M.setRegister(I.RD, M.getRegister(I.RS) - M.getRegister(I.RT));
    );

  IMPLEMENT(mul,
      int64_t Result;
      int32_t Half;
      Result =  (int32_t) M.getRegister(I.RS);
      Result *= (int32_t) M.getRegister(I.RT);

      Half = (Result & 0xFFFFFFFF);
      if (I.RD != 0)
        M.setRegister(I.RD, Half);

      Half = ((Result >> 32) & 0xFFFFFFFF);
      if (I.RV != 0)
        M.setRegister(I.RV, Half);
    );

  IMPLEMENT(mulu,
      uint64_t Result;
      int32_t Half;
      Result =  (uint32_t) M.getRegister(I.RS);
      Result *= (uint32_t) M.getRegister(I.RT);

      Half = (Result & 0xFFFFFFFF);
      if (I.RD != 0)
        M.setRegister(I.RD, Half);

      Half = ((Result >> 32) & 0xFFFFFFFF);
      if (I.RV != 0)
        M.setRegister(I.RV, Half);
    );

  IMPLEMENT(ext,
     uint32_t Lsb = I.RS;
     uint32_t Size = I.RT+1;
     M.setRegister(I.RV, (((unsigned)M.getRegister(I.RD)) << (32 - Size - Lsb)) >> (32 - Size));
   );

  IMPLEMENT(div,
      M.setRegister(I.RD, M.getRegister(I.RS) / M.getRegister(I.RT));
    );

  IMPLEMENT(mod,
      M.setRegister(I.RV, M.getRegister(I.RS) % M.getRegister(I.RT));
    );

  IMPLEMENT(divu,
      M.setRegister(I.RD, (uint32_t) M.getRegister(I.RS) / (uint32_t) M.getRegister(I.RT));
    );

  IMPLEMENT(modu,
      M.setRegister(I.RV, (uint32_t) M.getRegister(I.RS) % (uint32_t) M.getRegister(I.RT));
    );

  IMPLEMENT(ldhu,
      unsigned short int half = M.getMemHalfAt(M.getRegister(I.RS) + I.Imm);
      M.setRegister(I.RT, (uint32_t) half);
    );

  IMPLEMENT(ldh,
      short int half = M.getMemHalfAt(M.getRegister(I.RS) + I.Imm);
			M.setRegister(I.RT, (int32_t) half);
    );

  IMPLEMENT(sth,
      uint16_t half = M.getRegister(I.RT) & 0xFFFF;
      M.setMemByteAt(M.getRegister(I.RS) + I.Imm    , (half)      & 0xFF);
      M.setMemByteAt(M.getRegister(I.RS) + I.Imm + 1, (half >> 8) & 0xFF);
    );

  IMPLEMENT(ldi,
      M.setRegister(LDI_REG, I.RT);
      M.setRegister(I.RT, (M.getRegister(I.RT) & 0xFFFFC000) | (I.Imm & 0x3FFF));
    );

  IMPLEMENT(ldihi,
      M.setRegister(M.getRegister(LDI_REG), (M.getRegister(M.getRegister(LDI_REG)) & 0x3FFF) | (I.Addrs << 14));
    );

  IMPLEMENT(ldw,
      M.setRegister(I.RT, M.getMemValueAt(M.getRegister(I.RS) + I.Imm).asI_);
    );

  IMPLEMENT(addi,
      M.setRegister(I.RT, M.getRegister(I.RS) + I.Imm);
    );

  IMPLEMENT(and_,
      M.setRegister(I.RD, M.getRegister(I.RS) & M.getRegister(I.RT));
    );

  IMPLEMENT(andi,
      M.setRegister(I.RT, M.getRegister(I.RS) & (I.Imm & 0x3FFF));
    );

  IMPLEMENT(or_,
      M.setRegister(I.RD, M.getRegister(I.RS) | M.getRegister(I.RT));
    );

  IMPLEMENT(nor,
      M.setRegister(I.RD, ~(M.getRegister(I.RS) | M.getRegister(I.RT)));
    );

  IMPLEMENT(shr,
      unsigned aux = ((uint32_t) M.getRegister(I.RT)) >> (uint32_t) I.RS;
      M.setRegister(I.RD, aux);
    );

  IMPLEMENT(asr,
      M.setRegister(I.RD, ((int32_t) M.getRegister(I.RT)) >> I.RS);
    );

  IMPLEMENT(asrr,
      M.setRegister(I.RD, ((int32_t) M.getRegister(I.RT)) >> (M.getRegister(I.RS) & 0x1F));
    );

  IMPLEMENT(shl,
      M.setRegister(I.RD, M.getRegister(I.RT) << I.RS);
    );

  IMPLEMENT(shlr,
      M.setRegister(I.RD, M.getRegister(I.RT) << (M.getRegister(I.RS) & 0x1F));
    );

  IMPLEMENT(shrr,
      M.setRegister(I.RD, (uint32_t) M.getRegister(I.RT) >> (M.getRegister(I.RS) & 0x1F));
    );

  IMPLEMENT(movn,
      if (M.getRegister(I.RT) != 0)
        M.setRegister(I.RD, M.getRegister(I.RS));
    );

  IMPLEMENT(movz,
      if (M.getRegister(I.RT) == 0)
        M.setRegister(I.RD, M.getRegister(I.RS));
    );

  IMPLEMENT(ldb,
      M.setRegister(I.RT, (int32_t) M.getMemByteAt(M.getRegister(I.RS) + I.Imm));
    );

  IMPLEMENT(ldbu,
      M.setRegister(I.RT, (uint32_t) M.getMemByteAt(M.getRegister(I.RS) + I.Imm));
    );

  IMPLEMENT(seh,
      M.setRegister(I.RS, (int32_t) ((int16_t) M.getRegister(I.RT)));
    );

  IMPLEMENT(seb,
      M.setRegister(I.RS, (int32_t) ((int8_t) M.getRegister(I.RT)));
    );

  IMPLEMENT(stb,
      M.setMemByteAt(M.getRegister(I.RS) + I.Imm, (unsigned char) M.getRegister(I.RT) & 0xFF);
    );

  IMPLEMENT(stw,
      M.setMemValueAt(M.getRegister(I.RS) + I.Imm, M.getRegister(I.RT));
    );

  IMPLEMENT(sltiu,
      if (((uint32_t) M.getRegister(I.RS)) < ((uint32_t) (I.Imm & 0x3FFF)))
        M.setRegister(I.RT, 1);
      else
        M.setRegister(I.RT, 0);
    );

  IMPLEMENT(slti,
      M.setRegister(I.RT, (int32_t) M.getRegister(I.RS) < (int32_t) I.Imm);
    );

  IMPLEMENT(sltu,
      M.setRegister(I.RD, (uint32_t) M.getRegister(I.RS) < (uint32_t) M.getRegister(I.RT));
    );

  IMPLEMENT(slt,
      M.setRegister(I.RD, (int32_t) M.getRegister(I.RS) < (int32_t) M.getRegister(I.RT));
    );

  IMPLEMENT(xori,
      M.setRegister(I.RT, M.getRegister(I.RS) ^ (I.Imm & 0x3FFF));
    );

  IMPLEMENT(xor_,
      M.setRegister(I.RD, M.getRegister(I.RS) ^ M.getRegister(I.RT));
    );

  IMPLEMENT(ori,
      M.setRegister(I.RT, M.getRegister(I.RS) | (I.Imm & 0x3FFF));
    );

  IMPLEMENT(ror,
      M.setRegister(I.RD, rotate_right(M.getRegister(I.RT), I.RS));
    );

  IMPLEMENT(ijmphi,
      M.setRegister(IJMP_REG, 0 | (I.Addrs << 12));
    );

  /**********************  Float Inst  **************************/

   IMPLEMENT(absd,
       M.setDoubleRegister(I.RS, fabs(M.getDoubleRegister(I.RT)));
    );

   IMPLEMENT(abss,
       M.setFloatRegister(I.RS, fabs(M.getFloatRegister(I.RT)));
    );

   IMPLEMENT(ldc1,
			int32_t* R = M.getRegisterPtr();
			char* M1 = M.getByteMemoryPtr();	
			uint32_t Offset = M.getDataMemOffset();
			
			*(((int64_t*) R) + 65 + I.RT) = *((int64_t*) ((M1 + I.Imm + (*(R + I.RS))) - Offset));
    );

   IMPLEMENT(lwc1,
       M.setFloatRegister(I.RT, M.getMemValueAt(M.getRegister(I.RS) + I.Imm).asF_);
    );

   IMPLEMENT(lwxc1,
       M.setFloatRegister(I.RD, M.getMemValueAt(M.getRegister(I.RT) + M.getRegister(I.RS)).asF_);
    );

   IMPLEMENT(ldxc1,
			int32_t* R = M.getRegisterPtr();
			char* M1 = M.getByteMemoryPtr();	
			uint32_t Offset = M.getDataMemOffset();
			
			*(((int64_t*) R) + 65 + I.RD) = *((int64_t*) ((M1 + (*(R + I.RT)) + (*(R + I.RS))) - Offset));
    );

   IMPLEMENT(sdxc1,
      M.setMemValueAt(M.getRegister(I.RT) + M.getRegister(I.RS) +0, M.getRegister(130 + I.RD*2 + 0));
      M.setMemValueAt(M.getRegister(I.RT) + M.getRegister(I.RS) +4 , M.getRegister(130 + I.RD*2 + 1));
    );

   IMPLEMENT(sdc1,
      M.setMemValueAt(M.getRegister(I.RS) + I.Imm +0, M.getRegister(130 + I.RT*2+0));
      M.setMemValueAt(M.getRegister(I.RS) + I.Imm +4, M.getRegister(130 + I.RT*2+1));
    );

   IMPLEMENT(swc1,
      M.setMemValueAt(M.getRegister(I.RS) + I.Imm, M.getRegister(66 + I.RT));
    );

   IMPLEMENT(swxc1,
      M.setMemValueAt(M.getRegister(I.RT) + M.getRegister(I.RS), M.getRegister(66 + I.RD));
    );

   IMPLEMENT(mtlc1,
       DWordBit DW;
       DW.asF = M.getDoubleRegister(I.RT);
       DW.asI = (DW.asI & 0xFFFFFFFF00000000ULL) + (((uint64_t) M.getRegister(I.RS)));
       M.setDoubleRegister(I.RT, DW.asF);
    );

   IMPLEMENT(mthc1,
       DWordBit DW;
       DW.asF = M.getDoubleRegister(I.RT);
       DW.asI = (DW.asI & 0xFFFFFFFFULL) + (((uint64_t) M.getRegister(I.RS)) << 32);
       M.setDoubleRegister(I.RT, DW.asF);
    );

   IMPLEMENT(mflc1,
       DWordBit DW;
       DW.asF = M.getDoubleRegister(I.RT);
       M.setRegister(I.RS, (uint32_t)(DW.asI & 0xFFFFFFFF));
   );

   IMPLEMENT(mfhc1,
       DWordBit DW;
       DW.asF = M.getDoubleRegister(I.RT);
       M.setRegister(I.RS, (uint32_t)(DW.asI >> 32));
   );

   IMPLEMENT(ceqd,
       double A = M.getDoubleRegister(I.RS);
       double B = M.getDoubleRegister(I.RT);
       M.setRegister(CC_REG, A == B ? (isnan(A) || isnan(B) ? 0 : 1) : 0);
    );

   IMPLEMENT(ceqs,
       float A = M.getFloatRegister(I.RS);
       float B = M.getFloatRegister(I.RT);
       M.setRegister(CC_REG, A == B ? (isnan(A) || isnan(B) ? 0 : 1) : 0);
    );

   IMPLEMENT(negd,
       M.setDoubleRegister(I.RS, -M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(negs,
       M.setFloatRegister(I.RS, -M.getFloatRegister(I.RT));
    );

   IMPLEMENT(movd,
       M.setDoubleRegister(I.RS, M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(movf,
      if (M.getRegister(CC_REG) == 0)
        M.setRegister(I.RS, M.getRegister(I.RT));
    );

   IMPLEMENT(movt,
      if (M.getRegister(CC_REG) != 0)
        M.setRegister(I.RS, M.getRegister(I.RT));
    );

   IMPLEMENT(movts,
      if (M.getRegister(CC_REG) != 0)
        M.setFloatRegister(I.RS, M.getFloatRegister(I.RT));
    );
   IMPLEMENT(movs,
       M.setFloatRegister(I.RS, M.getFloatRegister(I.RT));
    );

   IMPLEMENT(movzd,
       if (M.getRegister(I.RT) == 0)
        M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS));
    );

   IMPLEMENT(movzs,
       if (M.getRegister(I.RT) == 0)
        M.setFloatRegister(I.RD, M.getFloatRegister(I.RS));
    );

   IMPLEMENT(movnd,
       if (M.getRegister(I.RT) != 0)
        M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS));
    );

   IMPLEMENT(movns,
       if (M.getRegister(I.RT) != 0)
        M.setFloatRegister(I.RD, M.getFloatRegister(I.RS));
    );

   IMPLEMENT(movfd,
       if (M.getRegister(CC_REG) == 0)
        M.setDoubleRegister(I.RS, M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(movfs,
       if (M.getRegister(CC_REG) == 0)
        M.setFloatRegister(I.RS, M.getFloatRegister(I.RT));
    );

   IMPLEMENT(movtd,
       if (M.getRegister(CC_REG) != 0)
        M.setDoubleRegister(I.RS, M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(adds,
       M.setFloatRegister(I.RD, M.getFloatRegister(I.RS) + M.getFloatRegister(I.RT));
    );

   IMPLEMENT(subd,
       M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS) - M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(subs,
       M.setFloatRegister(I.RD, M.getFloatRegister(I.RS) - M.getFloatRegister(I.RT));
    );

   IMPLEMENT(muls,
       M.setFloatRegister(I.RD, M.getFloatRegister(I.RS) * M.getFloatRegister(I.RT));
    );

   IMPLEMENT(muld,
       M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS) * M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(divd,
       M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS) / M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(divs,
       M.setFloatRegister(I.RD, M.getFloatRegister(I.RS) / M.getFloatRegister(I.RT));
    );

   IMPLEMENT(addd,
       M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS) + M.getDoubleRegister(I.RT));
   );

   IMPLEMENT(maddd,
       M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS) * M.getDoubleRegister(I.RT) + M.getDoubleRegister(I.RV));
   );

   IMPLEMENT(msubs,
       M.setFloatRegister(I.RD, M.getFloatRegister(I.RS) * M.getFloatRegister(I.RT) - M.getFloatRegister(I.RV));
   );

   IMPLEMENT(msubd,
       M.setDoubleRegister(I.RD, M.getDoubleRegister(I.RS) * M.getDoubleRegister(I.RT) - M.getDoubleRegister(I.RV));
   );

   IMPLEMENT(madds,
       M.setFloatRegister(I.RD, M.getFloatRegister(I.RS) * M.getFloatRegister(I.RT) + M.getFloatRegister(I.RV));
   );

   IMPLEMENT(mtc1,
       WordBit Tmp;
       Tmp.asI = M.getRegister(I.RS);
       M.setFloatRegister(I.RT, Tmp.asF);
    );

   IMPLEMENT(mfc1,
       WordBit Tmp;
       Tmp.asF = M.getFloatRegister(I.RT);
       M.setRegister(I.RS, Tmp.asI);
    );

   IMPLEMENT(truncws,
       WordBit Tmp;
       Tmp.asI = (int32_t) M.getFloatRegister(I.RT);
       M.setFloatRegister(I.RS, Tmp.asF);
    );

   IMPLEMENT(truncwd,
       WordBit Tmp;
       Tmp.asI = (int32_t) M.getDoubleRegister(I.RT);
       M.setFloatRegister(I.RS, Tmp.asF);
    );

   IMPLEMENT(cvtsw,
       WordBit Tmp;
       Tmp.asF = M.getFloatRegister(I.RT);
       M.setFloatRegister(I.RS, (float) (int) Tmp.asI);
    );

   IMPLEMENT(cvtdw,
       WordBit Tmp;
       Tmp.asF = M.getFloatRegister(I.RT);
       M.setDoubleRegister(I.RS, (double) (int) Tmp.asI);
    );

   IMPLEMENT(cvtds,
       M.setDoubleRegister(I.RS, (double) M.getFloatRegister(I.RT));
    );

   IMPLEMENT(cvtsd,
       M.setFloatRegister(I.RS, (float) M.getDoubleRegister(I.RT));
    );

   IMPLEMENT(coltd,
       double A = M.getDoubleRegister(I.RS);
       double B = M.getDoubleRegister(I.RT);
       M.setRegister(CC_REG, A < B ? (isnan(A) || isnan(B) ? 0 : 1) : 0);
    );

   IMPLEMENT(colts,
       double A = M.getFloatRegister(I.RS);
       double B = M.getFloatRegister(I.RT);
       M.setRegister(CC_REG, A < B ? (isnan(A) || isnan(B) ? 0 : 1) : 0);
    );

   IMPLEMENT(coled,
       double A = M.getDoubleRegister(I.RS);
       double B = M.getDoubleRegister(I.RT);
       M.setRegister(CC_REG, A <= B ? (isnan(A) || isnan(B) ? 0 : 1) : 0);
    );

   IMPLEMENT(coles,
       float A = M.getFloatRegister(I.RS);
       float B = M.getFloatRegister(I.RT);
       M.setRegister(CC_REG, A <= B ? (isnan(A) || isnan(B) ? 0 : 1) : 0);
    );

   IMPLEMENT(culed,
       M.setRegister(CC_REG, M.getDoubleRegister(I.RS) <= M.getDoubleRegister(I.RT) ? 1 : 0);
    );

   IMPLEMENT(cules,
       M.setRegister(CC_REG, M.getFloatRegister(I.RS) <= M.getFloatRegister(I.RT) ? 1 : 0);
    );

   IMPLEMENT(cults,
       M.setRegister(CC_REG, M.getFloatRegister(I.RS) < M.getFloatRegister(I.RT) ? 1 : 0);
    );

   IMPLEMENT(cultd,
       M.setRegister(CC_REG, M.getDoubleRegister(I.RS) < M.getDoubleRegister(I.RT) ? 1 : 0);
    );

   IMPLEMENT(cund,
       M.setRegister(CC_REG, (isnan(M.getDoubleRegister(I.RS)) || isnan(M.getDoubleRegister(I.RT))) ? 1 : 0);
    );

   IMPLEMENT(cuns,
       M.setRegister(CC_REG, (isnan(M.getFloatRegister(I.RS)) || isnan(M.getFloatRegister(I.RT))) ? 1 : 0);
    );

   IMPLEMENT(cueqs,
       M.setRegister(CC_REG, (M.getFloatRegister(I.RS) == M.getFloatRegister(I.RT)) ? 1 : 0);
    );

   IMPLEMENT(cueqd,
       M.setRegister(CC_REG, (M.getDoubleRegister(I.RS) == M.getDoubleRegister(I.RT)) ? 1 : 0);
    );

   IMPLEMENT(sqrtd,
       M.setDoubleRegister(I.RS, sqrt(M.getDoubleRegister(I.RT)));
    );

   IMPLEMENT(sqrts,
       M.setFloatRegister(I.RS, sqrt(M.getFloatRegister(I.RT)));
    );

  /********************** JMPs and BRs **************************/

  IMPLEMENT_BR(jeq,
      if (M.getRegister(I.RS) == M.getRegister(I.RT)) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jeqz,
      if (M.getRegister(I.RS) == 0) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jgtz,
      if (!(M.getRegister(I.RT) & 0x80000000) && (M.getRegister(I.RT) != 0)) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jgez,
      if (!(M.getRegister(I.RT) & 0x80000000)) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jlez,
      if ((M.getRegister(I.RT) == 0) || (M.getRegister(I.RT) & 0x80000000)) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jltz,
      if (M.getRegister(I.RT) & 0x80000000) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jne,
      if (M.getRegister(I.RS) != M.getRegister(I.RT)) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

  IMPLEMENT_BR(jnez,
      if (M.getRegister(I.RS) != 0) {
        M.setPC(M.getPC() + (I.Imm << 2) + 4);
        GOTO_NEXT;
      }
    );

   IMPLEMENT_BR(bc1f,
       if (M.getRegister(CC_REG) == 0) {
         M.setPC(M.getPC() + (I.Imm << 2) + 4);
         GOTO_NEXT;
       }
    );

   IMPLEMENT_BR(bc1t,
       if (M.getRegister(CC_REG) == 1) {
         M.setPC(M.getPC() + (I.Imm << 2) + 4);
         GOTO_NEXT;
       }
    );

  IMPLEMENT_JMP(call,
      M.setRegister(31, M.getPC()+4);
      M.setPC((M.getPC() & 0xF0000000) | (I.Addrs << 2));
      //#ifdef DEBUG
      //std::cerr << "Call to " << std::hex << M.getPC() <<
      //#endif
    );

  IMPLEMENT_JMP(callr,
      M.setRegister(31, M.getPC()+4);
      M.setPC(M.getRegister(I.RT));
    );

  IMPLEMENT_JMP(jumpr,
      M.setPC(M.getRegister(I.RT));
    );

  IMPLEMENT_JMP(jump,
      M.setPC((M.getPC() & 0xF0000000) | (I.Addrs << 2));
    );

  IMPLEMENT_JMP(ijmp,
      M.setRegister(IJMP_REG, M.getRegister(IJMP_REG) & 0xFFFFF000);
      M.setRegister(IJMP_REG, M.getRegister(IJMP_REG) | (I.Imm & 0xFFF));
      uint32_t Target = M.getMemValueAt(M.getRegister(IJMP_REG) + M.getRegister(I.RT)).asI_;
      M.setPC(Target);
    );

	IMPLEMENT(syscall,
    	if (SyscallM.processSyscall(M))
      	return;
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


