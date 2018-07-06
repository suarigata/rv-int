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
    case RVInstType::SLLI:
    case RVInstType::SRLI:
    case RVInstType::SRAI:
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
    case RVInstType::FADD.S:
    case RVInstType::FSUB.S:
    case RVInstType::FMUL.S:
    case RVInstType::FDIV.S:
    case RVInstType::FSQRT.S:
    case RVInstType::FSGNJ.S:
    case RVInstType::FSGNJN.S:
    case RVInstType::FSGNJX.S:
    case RVInstType::FMIN.S:
    case RVInstType::FMAX.S:
    case RVInstType::FCVT.W.S:
    case RVInstType::FCVT.WU.S:
    case RVInstType::FMV.X.W:
    case RVInstType::FEQ.S:
    case RVInstType::FLT.S:
    case RVInstType::FLE.S:
    case RVInstType::FCLASS.S:
    case RVInstType::FCVT.S.W:
    case RVInstType::FCVT.S.WU:
    case RVInstType::FMV.W.X:
    case RVInstType::FADD.D:
    case RVInstType::FSUB.D:
    case RVInstType::FMUL.D:
    case RVInstType::FDIV.D:
    case RVInstType::FSQRT.D:
    case RVInstType::FSGNJ.D:
    case RVInstType::FSGNJN.D:
    case RVInstType::FSGNJX.D:
    case RVInstType::FMIN.D:
    case RVInstType::FMAX.D:
    case RVInstType::FCVT.S.D:
    case RVInstType::FCVT.D.S:
    case RVInstType::FEQ.D:
    case RVInstType::FLT.D:
    case RVInstType::FLE.D:
    case RVInstType::FCLASS.D:
    case RVInstType::FCVT.W.D:
    case RVInstType::FCVT.WU.D:
    case RVInstType::FCVT.D.W:
    case RVInstType::FCVT.D.WU:
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
    case RVInstType::FENCE:
    case RVInstType::FENCE.I:
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
    case RVInstType::FMADD.S:
    case RVInstType::FMSUB.S:
    case RVInstType::FNMSUB.S:
    case RVInstType::FNMADD.S:
    case RVInstType::FMADD.D:
    case RVInstType::FMSUB.D:
    case RVInstType::FNMSUB.D:
    case RVInstType::FNMADD.D:
      return EncodingType::ER4;
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
          if(getF7(W))
            I.Type = RVInstType::SRAI;
          else
            I.Type = RVInstType::SRLI;
          break;
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
            if(getF7(W))
              I.Type = RVInstType::SUB;
            else
              I.Type = RVInstType::ADD;
            break;
          case 0b001: I.Type = RVInstType::SLL;		break; // SLL
          case 0b010: I.Type = RVInstType::SLT;		break; // SLT
          case 0b011: I.Type = RVInstType::SLTU;	break; // SLTU
          case 0b100: I.Type = RVInstType::XOR;		break; // XOR
          case 0b101: // SRL | SRA
            if(getF7(W))
              I.Type = RVInstType::SRA;
            else
              I.Type = RVInstType::SRL;
            break;
          case 0b110: I.Type = RVInstType::OR;		break; // OR
          case 0b111: I.Type = RVInstType::AND;		break; // AND
        }
      break;
    case 0b0001111: // FENCE | FENCE.I
      I.Type = RVInstType::ADD; // TODO
      break;
    case 0b1110011: // SYSTEM
      switch(getF3(F)){
        case 0b000: // ECALL | EBREAK
          if(getImmI(W))
            I.Type = RVInstType::EBREAK;
          else
            I.Type = RVInstType::ECALL;
          break;
        case 0b001: I.Type = RVInstType::CSRRW;		break; // CSRRW
        case 0b010: I.Type = RVInstType::CSRRS;		break; // CSRRS
        case 0b011: I.Type = RVInstType::CSRRC;		break; // CSRRC
        case 0b101: I.Type = RVInstType::CSRRWI;	break; // CSRRWI
        case 0b110: I.Type = RVInstType::CSRRSI;	break; // CSRRSI
        case 0b111: I.Type = RVInstType::CSRRCI;	break; // CSRRCI
      }
      break;

    default:
      I.Type = RVInstType::Null;
  }

  if (I.Type == RVInstType::Null) {
    std::cout << "Houston: we have a problem! Inst (" << std::hex << CodedInst << ") not implemented! \n";
    exit(1);
  }

  dbt::RVDecoder::fillFields(I, dbt::RVDecoder::getEncodingType(I.Type), W);

  if (!CallZero && I.Type == Call && (I.Addrs << 2) == 0) { // TODO Call?
    std::cerr << "Pay attention! Something must have been linked wrongly, a call 0 was found! Maybe a -lm missing?\n";
    CallZero = true;
  }

  return I;
}
