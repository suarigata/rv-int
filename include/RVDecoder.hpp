#ifndef RVDECODER_HPP
#define RVDECODER_HPP

#include <memory>
#include <iostream>

#define LDI_REG   64
#define IJMP_REG  65
#define CC_REG    257

//#define PRINTREG

namespace dbt {
  namespace RVDecoder {
    union Word {
      char asC_[4];
      uint32_t asI_;
    };

    enum RVInstType{
      LUI, AUIPC, JAL, JALR, BEQ, BNE, BLT, BGE, BLTU, BGEU, LB, LH, LW, LBU, LHU, SB, SH, SW, ADDI, SLTI, SLTIU, XORI,
      ORI, ANDI, SLLI, SRLI, SRAI, ADD, SUB, SLL, SLT, SLTU, XOR, SRL, SRA, OR, AND, FENCE, FENCE.I, ECALL, EBREAK,
      CSRRW, CSRRS, CSRRC, CSRRWI, CSRRSI, CSRRCI,
      MUL, MULH, MULHSU, MULHU, DIV, DIVU, REM, REMU,
      FLW, FSW, FMADD.S, FMSUB.S, FNMSUB.S, FNMADD.S, FADD.S, FSUB.S, FMUL.S, FDIV.S, FSQRT.S, FSGNJ.S, FSGNJN.S, FSGNJX.S,
      FMIN.S, FMAX.S, FCVT.W.S, FCVT.WU.S, FMV.X.W, FEQ.S, FLT.S, FLE.S, FCLASS.S, FCVT.S.W, FCVT.S.WU, FMV.W.X,
      FLD, FSD, FMADD.D, FMSUB.D, FNMSUB.D, FNMADD.D, FADD.D, FSUB.D, FMUL.D, FDIV.D, FSQRT.D, FSGNJ.D, FSGNJN.D, FSGNJX.D,
      FMIN.D, FMAX.D, FCVT.S.D, FCVT.D.S, FEQ.D, FLT.D, FLE.D, FCLASS.D, FCVT.W.D, FCVT.WU.D, FCVT.D.W, FCVT.D.WU,
      Null
/*
      Add  , And , Andi, Or  , Ldi , Ldihi, Ldw , Addi, Call , Jumpr, Stw , Sltiu  , Slti   , Jeq   , Jne  , Jump , Mul  , 
      Shr  , Shl , Jeqz, Sub , Slt , Div  , Mod , Ori , Jgtz , Jlez , Jnez, Ldbu   , Stb    , Sltu  , Asr  , Jltz , Movn ,
      Nor  , Ldh , Ldb , Sth , Ldhu, Jgez , Nop , Seh , Callr, Shlr , Xor , Seb    , Ijmphi , Ijmp  , Divu , Modu , Ldc1 , 
      Mthc1, Ceqd, Ceqs, Bc1f, Bc1t, Movd , Lwc1, Adds, Addd , Mtc1 , Mfc1, Truncws, Truncwd, Cvtsw , Cvtdw, Cvtds, Cvtsd, 
      Mulu , Movz, Xori, Sdc1, Swc1, Maddd, Movs, Muls, Coltd, Swxc1, Negd, Lwxc1  , Syscall, Mtlc1 , Divs , Subs , Mflc1, 
      Mfhc1, Divd, Subd, Negs, Ext , Madds, Shrr, Movf, Movt , Ldxc1, Muld, Sdxc1  , Msubs  , Coled , Culed, Msubd, Movzd,
      Movfd, Asrr, Absd, Abss, Cund, Movnd, Ror, Movzs, Movfs, Colts, Movns, Coles , Sqrtd  , Sqrts , Cults, Cules, Cultd,
      Movtd, Movts,Cuns, Cueqs, Cueqd, Null
//*/
    };

    enum EncodingType {
      ER, ER4, EI, ES, EB, EU, EJ
//      PL0, PL6, PL26ij, PL26j, PL26c, PL26i, PL12, PL18, PL16, PL24, PL18i, PL20, PL20i
    };

    typedef struct RVInst {
      RVInstType Type;
      uint8_t RD, RS1, RS2, RS3, RM;
      int32_t Imm;
    } RVInst;

    EncodingType getEncodingType(RVInstType);
    void fillFields(RVInst&, EncodingType, Word); 
    RVInst decode(uint32_t);
    bool isControlFlowInst(RVInst);
    bool isIndirectBranch(RVInst);
    std::array<uint32_t, 2> getPossibleTargets(uint32_t, RVInst);
  }
}

#endif
