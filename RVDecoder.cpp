#include "RVDecoder.hpp"

using namespace dbt::RVDecoder;

constexpr uint32_t getImmU(Word W){
  return W.asI_ & 0xFFFFF000;
}

constexpr uint32_t getImmI(Word W){
  return W.asI_ >> 20;
}

constexpr uint32_t getImmS(Word W){
  return ((W.asI_ & 0xFE000000) >> 20) | ((W.asI_ >> 7) & 0x1F);
}

constexpr uint32_t getImmB(Word W){
  return ((W.asI_ & 0x80000000) >> 19) | ((W.asI_ & 0x80) << 4) | ((W.asI_ & 0x7E000000) >> 20) | ((W.asI_ >> 7) & 0x1E);
}

constexpr uint32_t getImmJ(Word W){
  return ((W.asI_ >> 11) & 0x100000) | (W.asI_ & 0xFF000) | ((W.asI_ >> 9) & 0x800) | ((W.asI_ >> 20) & 0x7FE);
}

constexpr uint8_t getRD(Word W){
  return (W.asI_ >> 7) & 0x1F;
}

constexpr uint8_t getRS1(Word W){
  return (W.asI_ >> 15) & 0x1F;
}

constexpr uint8_t getRS2(Word W){
  return (W.asI_ >> 20) & 0x1F;
}

constexpr uint8_t getRS3(Word W){
  return (W.asI_ >> 27); // & 0x1F;
}

constexpr uint8_t getF2(Word W){
  return (W.asI_ >> 25) & 0x3;
}

constexpr uint8_t getF3(Word W){
  return (W.asI_ >> 12) & 0x7;
}

constexpr uint8_t getF7(Word W){
  return (W.asI_ >> 25); // & 0x7F;
}

EncodingType dbt::RVDecoder::getEncodingType(RVInstType InstType) {
  switch(InstType) {
    case RVInstType::ADD:
    case RVInstType::SUB:
    case RVInstType::SLL:
    case RVInstType::SLT:
    case RVInstType::SLTU:
    case RVInstType::XOR:
    case RVInstType::SRL:
    case RVInstType::SRA:
    case RVInstType::OR:
    case RVInstType::AND:
    case RVInstType::MUL:
    case RVInstType::MULH:
    case RVInstType::MULHSU:
    case RVInstType::MULHU:
    case RVInstType::DIV:
    case RVInstType::DIVU:
    case RVInstType::REM:
    case RVInstType::REMU:
    case RVInstType::FADD_S:
    case RVInstType::FSUB_S:
    case RVInstType::FMUL_S:
    case RVInstType::FDIV_S:
    case RVInstType::FSQRT_S:
    case RVInstType::FSGNJ_S:
    case RVInstType::FSGNJN_S:
    case RVInstType::FSGNJX_S:
    case RVInstType::FMIN_S:
    case RVInstType::FMAX_S:
    case RVInstType::FCVT_W_S:
    case RVInstType::FCVT_WU_S:
    case RVInstType::FMV_X_W:
    case RVInstType::FEQ_S:
    case RVInstType::FLT_S:
    case RVInstType::FLE_S:
    case RVInstType::FCLASS_S:
    case RVInstType::FCVT_S_W:
    case RVInstType::FCVT_S_WU:
    case RVInstType::FMV_W_X:
    case RVInstType::FADD_D:
    case RVInstType::FSUB_D:
    case RVInstType::FMUL_D:
    case RVInstType::FDIV_D:
    case RVInstType::FSQRT_D:
    case RVInstType::FSGNJ_D:
    case RVInstType::FSGNJN_D:
    case RVInstType::FSGNJX_D:
    case RVInstType::FMIN_D:
    case RVInstType::FMAX_D:
    case RVInstType::FCVT_S_D:
    case RVInstType::FCVT_D_S:
    case RVInstType::FEQ_D:
    case RVInstType::FLT_D:
    case RVInstType::FLE_D:
    case RVInstType::FCLASS_D:
    case RVInstType::FCVT_W_D:
    case RVInstType::FCVT_WU_D:
    case RVInstType::FCVT_D_W:
    case RVInstType::FCVT_D_WU:
      return EncodingType::ER;
    case RVInstType::JALR:
    case RVInstType::LB:
    case RVInstType::LH:
    case RVInstType::LW:
    case RVInstType::LBU:
    case RVInstType::LHU:
    case RVInstType::ADDI:
    case RVInstType::SLTI:
    case RVInstType::SLTIU:
    case RVInstType::XORI:
    case RVInstType::ORI:
    case RVInstType::ANDI:
    case RVInstType::SLLI:
    case RVInstType::SRLI:
    case RVInstType::SRAI:
    case RVInstType::FENCE:
    case RVInstType::FENCE_I:
    case RVInstType::ECALL:
    case RVInstType::EBREAK:
    case RVInstType::CSRRW:
    case RVInstType::CSRRS:
    case RVInstType::CSRRC:
    case RVInstType::CSRRWI:
    case RVInstType::CSRRSI:
    case RVInstType::CSRRCI:
    case RVInstType::FLW:
    case RVInstType::FLD:
      return EncodingType::EI;
    case RVInstType::SB:
    case RVInstType::SH:
    case RVInstType::SW:
    case RVInstType::FSW:
    case RVInstType::FSD:
      return EncodingType::ES;
    case RVInstType::BEQ:
    case RVInstType::BNE:
    case RVInstType::BLT:
    case RVInstType::BGE:
    case RVInstType::BLTU:
    case RVInstType::BGEU:
      return EncodingType::EB;
    case RVInstType::LUI:
    case RVInstType::AUIPC:
      return EncodingType::EU;
    case RVInstType::JAL:
      return EncodingType::EJ;
    case RVInstType::FMADD_S:
    case RVInstType::FMSUB_S:
    case RVInstType::FNMSUB_S:
    case RVInstType::FNMADD_S:
    case RVInstType::FMADD_D:
    case RVInstType::FMSUB_D:
    case RVInstType::FNMSUB_D:
    case RVInstType::FNMADD_D:
      return EncodingType::ER4;
    default:;
  }

  std::cout << "Dammit! We have a bug on Encoding types!\n";
  return EncodingType::ER;
}

extern void dbt::RVDecoder::fillFields(RVInst& I, EncodingType E, Word W) {
  switch (E) {
    case EncodingType::ER:
      I.RD = getRD(W);
      I.RS1 = getRS1(W);
      I.RS2 = getRS2(W);
      I.RM = getF3(W);
      break;
    case EncodingType::EI:
      I.Imm = getImmI(W);
      I.RD = getRD(W);
      I.RS1 = getRS1(W);
      break;
    case EncodingType::ES:
      I.Imm = getImmS(W);
      I.RS1 = getRS1(W);
      I.RS2 = getRS2(W);
      break;
    case EncodingType::EB:
      I.Imm = getImmB(W);
      I.RS1 = getRS1(W);
      I.RS2 = getRS2(W);
      break;
    case EncodingType::EU:
      I.Imm = getImmU(W);
      I.RD = getRD(W);
      break;
    case EncodingType::EJ:
      I.Imm = getImmJ(W);
      I.RD = getRD(W);
      break;
    case EncodingType::ER4:
      I.RD = getRD(W);
      I.RS1 = getRS1(W);
      I.RS2 = getRS2(W);
      I.RS3 = getRS3(W);
      I.RM = getF3(W);
      break;

#ifdef PRINTREG
      //std::cout << "RS: " << I.RS << "; RT: " << I.RT << "; Imm: " << I.Imm << ";\n";
#endif

    default: 
      std::cout << "It may seem unbelievable, but we haven't implemented this encoding type!\n";
      break;
  }
}

bool CallZero;
RVInst dbt::RVDecoder::decode(uint32_t CodedInst) {
  Word W;
  W.asI_ = CodedInst;
  uint8_t Op = W.asI_ & 0x7F;
  RVInst I;

//I.RD=I.RS1=I.RS2=I.RS3=I.RM=I.Imm=0; // TODO pode tirar



  I.Type = RVInstType::Null;
  
  switch(Op){
    case 0b0110111: I.Type = RVInstType::LUI;	break; // LUI
    case 0b0010111: I.Type = RVInstType::AUIPC;	break; // AUIPC
    case 0b1101111: I.Type = RVInstType::JAL;	break; // JAL
    case 0b1100111: I.Type = RVInstType::JALR;	break; // JALR
    case 0b1100011: // Branches
      switch(getF3(W)){
        case 0b000: I.Type = RVInstType::BEQ;	break; // BEQ
        case 0b001: I.Type = RVInstType::BNE;	break; // BNE
        case 0b100: I.Type = RVInstType::BLT;	break; // BLT
        case 0b101: I.Type = RVInstType::BGE;	break; // BGE
        case 0b110: I.Type = RVInstType::BLTU;	break; // BLTU
        case 0b111: I.Type = RVInstType::BGEU;	break; // BGEU
      }
      break;
    case 0b0000011: // Loads
      switch(getF3(W)){
        case 0b000: I.Type = RVInstType::LB;	break; // LB
        case 0b001: I.Type = RVInstType::LH;	break; // LH
        case 0b010: I.Type = RVInstType::LW;	break; // LW
        case 0b100: I.Type = RVInstType::LBU;	break; // LBU
        case 0b101: I.Type = RVInstType::LHU;	break; // LHU
      }
      break;
    case 0b0100011: // Stores
      switch(getF3(W)){
        case 0b000: I.Type = RVInstType::SB;	break; // SB
        case 0b001: I.Type = RVInstType::SH;	break; // SH
        case 0b010: I.Type = RVInstType::SW;	break; // SW
      }
      break;
    case 0b0010011:
      switch(getF3(W)){
        case 0b000: I.Type = RVInstType::ADDI;	break; // ADDI
        case 0b010: I.Type = RVInstType::SLTI;	break; // SLTI
        case 0b011: I.Type = RVInstType::SLTIU;	break; // SLTIU
        case 0b100: I.Type = RVInstType::XORI;	break; // XORI
        case 0b110: I.Type = RVInstType::ORI;	break; // ORI
        case 0b111: I.Type = RVInstType::ANDI;	break; // ANDI
        case 0b001: I.Type = RVInstType::SLLI;	break; // SLLI
        case 0b101: // SRLI | SRAI
          I.Type = (getF7(W) ? RVInstType::SRAI : RVInstType::SRLI); break;
      }
      break;
    case 0b0110011:
      if(getF7(W) & 0x1) // M Extension
        switch(getF3(W)){
          case 0b000: I.Type = RVInstType::MUL;		break; // MUL
          case 0b001: I.Type = RVInstType::MULH;	break; // MULH
          case 0b010: I.Type = RVInstType::MULHSU;	break; // MULHSU
          case 0b011: I.Type = RVInstType::MULHU;	break; // MULHU
          case 0b100: I.Type = RVInstType::DIV;		break; // DIV
          case 0b101: I.Type = RVInstType::DIVU;	break; // DIVU
          case 0b110: I.Type = RVInstType::REM;		break; // REM
          case 0b111: I.Type = RVInstType::REMU;	break; // REMU
        }
      else
        switch(getF3(W)){
          case 0b000: // ADD | SUB
            I.Type = (getF7(W) ? RVInstType::SUB : RVInstType::ADD); break;
          case 0b001: I.Type = RVInstType::SLL;		break; // SLL
          case 0b010: I.Type = RVInstType::SLT;		break; // SLT
          case 0b011: I.Type = RVInstType::SLTU;	break; // SLTU
          case 0b100: I.Type = RVInstType::XOR;		break; // XOR
          case 0b101: // SRL | SRA
            I.Type = (getF7(W) ? RVInstType::SRA : RVInstType::SRL); break;
          case 0b110: I.Type = RVInstType::OR;		break; // OR
          case 0b111: I.Type = RVInstType::AND;		break; // AND
        }
      break;
    case 0b0001111: // FENCE | FENCE.I
      I.Type = RVInstType::ADD; // TODO
      break;
    case 0b1110011: // SYSTEM
      switch(getF3(W)){
        case 0b000: // ECALL | EBREAK
          I.Type = ((getImmI(W)) ? RVInstType::EBREAK : RVInstType::ECALL); break;
        case 0b001: I.Type = RVInstType::CSRRW;		break; // CSRRW
        case 0b010: I.Type = RVInstType::CSRRS;		break; // CSRRS
        case 0b011: I.Type = RVInstType::CSRRC;		break; // CSRRC
        case 0b101: I.Type = RVInstType::CSRRWI;	break; // CSRRWI
        case 0b110: I.Type = RVInstType::CSRRSI;	break; // CSRRSI
        case 0b111: I.Type = RVInstType::CSRRCI;	break; // CSRRCI
      }
      break;
    case 0b0000111:
      I.Type = ((getF3(W) & 0x1) ? RVInstType::FLD : RVInstType::FLW); break;
    case 0b0100111:
      I.Type = ((getF3(W) & 0x1) ? RVInstType::FSD : RVInstType::FSW); break;
    case 0b1000011:
      I.Type = ((getF2(W)) ? RVInstType::FMADD_D : RVInstType::FMADD_S); break;
    case 0b1000111:
      I.Type = ((getF2(W)) ? RVInstType::FMSUB_D : RVInstType::FMSUB_S); break;
    case 0b1001011:
      I.Type = ((getF2(W)) ? RVInstType::FNMSUB_D : RVInstType::FNMSUB_S); break;
    case 0b1001111:
      I.Type = ((getF2(W)) ? RVInstType::FNMADD_D : RVInstType::FNMADD_S); break;
    case 0b1010011:
      switch(getF7(W)){
        case 0b0000000: I.Type = RVInstType::FADD_S;	break; // FADD_S
        case 0b0000001: I.Type = RVInstType::FADD_D;	break; // FADD_D
        case 0b0000100: I.Type = RVInstType::FSUB_S;	break; // FSUB_S
        case 0b0000101: I.Type = RVInstType::FSUB_D;	break; // FSUB_D
        case 0b0001000: I.Type = RVInstType::FMUL_S;	break; // FMUL_S
        case 0b0001001: I.Type = RVInstType::FMUL_D;	break; // FMUL_D
        case 0b0001100: I.Type = RVInstType::FDIV_S;	break; // FDIV_S
        case 0b0001101: I.Type = RVInstType::FDIV_D;	break; // FDIV_D
        case 0b0010000:
          switch(getF3(W)){
            case 0b000: I.Type = RVInstType::FSGNJ_S;	break; // FSGNJ_S
            case 0b001: I.Type = RVInstType::FSGNJN_S;	break; // FSGNJN_S
            case 0b010: I.Type = RVInstType::FSGNJX_S;	break; // FSGNJX_S
          }
          break;
        case 0b0010001:
          switch(getF3(W)){
            case 0b000: I.Type = RVInstType::FSGNJ_D;	break; // FSGNJ_D
            case 0b001: I.Type = RVInstType::FSGNJN_D;	break; // FSGNJN_D
            case 0b010: I.Type = RVInstType::FSGNJX_D;	break; // FSGNJX_D
          }
          break;
        case 0b0010100:
          I.Type = (getF3(W) ? RVInstType::FMAX_S : RVInstType::FMIN_S); break;
        case 0b0010101:
          I.Type = (getF3(W) ? RVInstType::FMAX_D : RVInstType::FMIN_D); break;
        case 0b0100000: I.Type = RVInstType::FCVT_S_D;	break; // FCVT_S_D
        case 0b0100001: I.Type = RVInstType::FCVT_D_S;	break; // FCVT_D_S
        case 0b0101100: I.Type = RVInstType::FSQRT_S;	break; // FSQRT_S
        case 0b0101101: I.Type = RVInstType::FSQRT_D;	break; // FSQRT_D
        case 0b1010000:
          switch(getF3(W)){
            case 0b000: I.Type = RVInstType::FLE_S;	break; // FLE_S
            case 0b001: I.Type = RVInstType::FLT_S;	break; // FLT_S
            case 0b010: I.Type = RVInstType::FEQ_S;	break; // FEQ_S
          }
          break;
        case 0b1010001:
          switch(getF3(W)){
            case 0b000: I.Type = RVInstType::FLE_D;	break; // FLE_D
            case 0b001: I.Type = RVInstType::FLT_D;	break; // FLT_D
            case 0b010: I.Type = RVInstType::FEQ_D;	break; // FEQ_D
          }
          break;
        case 0b1100000:
          I.Type = (getRS2(W) ? RVInstType::FCVT_WU_S : RVInstType::FCVT_W_S); break;
        case 0b1100001:
          I.Type = (getRS2(W) ? RVInstType::FCVT_WU_D : RVInstType::FCVT_W_D); break;
        case 0b1101000:
          I.Type = (getRS2(W) ? RVInstType::FCVT_S_WU : RVInstType::FCVT_S_W); break;
        case 0b1101001:
          I.Type = (getRS2(W) ? RVInstType::FCVT_D_WU : RVInstType::FCVT_D_W); break;
        case 0b1110000:
          I.Type = (getF3(W) ? RVInstType::FCLASS_S : RVInstType::FMV_X_W); break;
        case 0b1110001: I.Type = RVInstType::FCLASS_D;	break; // FCLASS_D
        case 0b1111000: I.Type = RVInstType::FMV_W_X;	break; // FMV_W_X
      }
      break;
  }

  if (I.Type == RVInstType::Null) {
    std::cout << "Houston: we have a problem! Inst (" << std::hex << CodedInst << ") not implemented!\n";
    exit(1);
  }
  
  dbt::RVDecoder::fillFields(I, dbt::RVDecoder::getEncodingType(I.Type), W);

/* TODO ver se isso é importante
  if (!CallZero && I.Type == Call && (I.Addrs << 2) == 0) { // TODO Call?
    std::cerr << "Pay attention! Something must have been linked wrongly, a call 0 was found! Maybe a -lm missing?\n";
    CallZero = true;
  }
//*/
  return I;
}

