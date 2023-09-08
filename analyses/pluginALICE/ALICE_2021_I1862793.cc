// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  /// @brief Prompt charm-strange baryons in pp at 13 TeV
  class ALICE_2021_I1862793 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1862793);

    void init() {

      const UnstableParticles up(Cuts::absrap < 0.5);
      declare(up, "up");

      book(_h_X0,1,1,1);                                     // XiC0 production cross section
      book(_h_XP,2,1,1);                                     // XiC+ production cross section
      book(_h_X0D0,3,1,1);                                   // ratio of production cross section between XiC0 and D0
      book(_h_XPDP,4,1,1);                                   // ratio of production cross section between XiC+ and D0
      book(_h_X0Lc,5,1,1);                                   // ratio of production cross section between XiC and LambdaC
      book(_h_XcSc,6,1,1);                                   // ratio of production cross section between XiC and SigmaC
      book(_h_X0Iint,7,1,1);                                 // XiC0 pt-integrated production cross section with 1<pt<12 GeV/c
      book(_h_X0Oint,8,1,1);                                 // XiC0 pt-integrated production cross section with pt>0 GeV/c
      book(_h_XPint,9,1,1);                                  // XiC+ pt-integrated production cross section with 4<pt<12 GeV/c


      book(_h_X0D,"TMP/_h_X0D",refData(3,1,1));              // XiC0 production cross section
      book(_h_D0,"TMP/_h_D0",refData(3,1,1));                // D0 production cross section

      book(_h_DP,"TMP/_h_DP",refData(4,1,1));                // D0 production cross section
      book(_h_XPD,"TMP/_h_XPD",refData(4,1,1));              // XiC+ production cross section

      book(_h_Lc,"TMP/_h_Lc",refData(5,1,1));                // LambdaC production cross section
      book(_h_X0L,"TMP/_h_X0L",refData(5,1,1));              // XiC0 production cross section

      book(_h_Sc,"TMP/_h_Sc",refData(6,1,1));                // SigmaC production cross section
      book(_h_Xc,"TMP/_h_Xc",refData(6,1,1));                // XiC production cross section


    }


    void analyze(const Event& event) {

      const UnstableParticles& up = apply<UnstableParticles>(event, "up");

      for (const Particle& p : up.particles()) {
        if(p.fromBottom())
          continue;

        else{
          if(p.abspid() == 4222 || p.abspid() == 4212 || p.abspid() == 4112){
            _h_Sc->fill(p.pT()/GeV);
          }

          else if(p.abspid() == 4122){
            _h_Lc->fill(p.pT()/GeV);
          }

          else if(p.abspid() == 421){
            _h_D0->fill(p.pT()/GeV);
            _h_DP->fill(p.pT()/GeV);
          }

          else if(p.abspid() == 4132){
            _h_X0->fill(p.pT()/GeV);
            _h_X0D->fill(p.pT()/GeV);
            _h_Xc->fill(p.pT()/GeV);
            _h_X0L->fill(p.pT()/GeV);
            _h_X0Oint->fill(1);

            if(p.pT()/GeV < 12.0 && p.pT()/GeV >= 1.0){
              _h_X0Iint->fill(1);
            }

          }

          else if(p.abspid() == 4232){
            _h_XP->fill(p.pT()/GeV);
            _h_Xc->fill(p.pT()/GeV);
            _h_XPD->fill(p.pT()/GeV);

             if(p.pT()/GeV < 12.0 && p.pT()/GeV >= 4.0){
              _h_XPint->fill(8);
            }
          }
        }
      }
    }

    void finalize() {

      scale(_h_Sc,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_Lc,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_D0,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_DP,              crossSection()/(microbarn*2*sumOfWeights()));

      scale(_h_Xc,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_X0,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_XP,              crossSection()/(microbarn*2*sumOfWeights()));

      scale(_h_X0D,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_XPD,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_X0L,              crossSection()/(microbarn*2*sumOfWeights()));

      divide(_h_X0D, _h_D0, _h_X0D0);
      divide(_h_XPD, _h_DP, _h_XPDP);
      divide(_h_X0L, _h_Lc, _h_X0Lc);
      divide(_h_Xc, _h_Sc, _h_XcSc);

      scale(_h_X0Iint,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_X0Oint,              crossSection()/(microbarn*2*sumOfWeights()));
      scale(_h_XPint,               crossSection()/(microbarn*2*sumOfWeights()));

    }

    Histo1DPtr _h_X0, _h_XP, _h_X0Iint,_h_X0Oint, _h_XPint, _h_Sc, _h_Lc, _h_D0, _h_DP, _h_Xc, _h_X0D, _h_XPD, _h_X0L;
    Scatter2DPtr _h_X0D0, _h_XPDP, _h_X0Lc, _h_XcSc;

  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1862793);

}
