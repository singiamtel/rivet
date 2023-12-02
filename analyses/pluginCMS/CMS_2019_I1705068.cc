// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// W + charm quark in proton-proton collisions at 13 TeV
  class CMS_2019_I1705068 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2019_I1705068);


    /// Initialise analyis
    void init() {

      // Projections
      declare(MissingMom(), "MET");
      LeptonFinder lf(0.1, Cuts::abseta < 2.4 && Cuts::abspid == PID::MUON);
      declare(lf, "Leptons");

      UnstableParticles dst(Cuts::pT > 5*GeV && Cuts::abseta < 2.4);
      declare(dst, "Dstar");

      // Particle-Level histograms from the paper:
      book(_hist_WplusMinus_MuAbseta, "d04-x01-y01");
      book(_hist_Wplus_MuAbseta, "d05-x01-y01");
      book(_hist_Wminus_MuAbseta, "d06-x01-y01");
    }


    void analyze(const Event& event) {

      // Identify the closest-matching mu+MET to m == mW
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      const Particles& mus = apply<LeptonFinder>(event, "Leptons").particles();
      const Particles mus_mtfilt = select(mus, [&](const Particle& m){ return mT(m, pmiss) > 0*GeV; });
      const int imfound = closestMatchIndex(mus_mtfilt, pmiss, Kin::mass, 80.4*GeV);

      // No MET or MT cut at generator level
      if (imfound < 0) vetoEvent;
      const Particle& lepton0 = mus_mtfilt[imfound];
      const double pt0 = lepton0.pT();
      const double eta0 = lepton0.abseta();
      const int muID = lepton0.pid();
      if (eta0 > 2.4 || pt0 < 26*GeV) vetoEvent;

      // D* selection:
      // OS = W boson and D* Meson have Opposite (charge) Signs
      // SS = W Boson and D* Meson have Same (charge) Signs
      // Associated W+c only has OS contributions,
      // W+ccbar (ccbar from gluon splitting) has equal probability to be OS or SS
      // OS-SS to remove the gluon splitting background

      const UnstableParticles& dst = apply<UnstableParticles>(event, "Dstar");
      for (const Particle& p : dst.particles()) {
        if (muID == -13 && p.pid() == -413) { // OS
          _hist_Wplus_MuAbseta->fill(eta0);
          _hist_WplusMinus_MuAbseta->fill(eta0);
        }
        else if (muID == 13 && p.pid() == 413) { // OS
          _hist_Wminus_MuAbseta->fill(eta0);
          _hist_WplusMinus_MuAbseta->fill(eta0);
        }
        else if (muID == -13 && p.pid() == 413) { // SS
          _hist_Wplus_MuAbseta->fill(eta0*-1);
          _hist_WplusMinus_MuAbseta->fill(eta0*-1);
        }
        else if (muID == 13 && p.pid() == -413) { // SS
          _hist_Wminus_MuAbseta->fill(eta0*-1);
          _hist_WplusMinus_MuAbseta->fill(eta0*-1);
        }
      }

    }


    void finalize() {
      scale(_hist_Wplus_MuAbseta, crossSection()/picobarn/sumOfWeights());
      scale(_hist_Wminus_MuAbseta, crossSection()/picobarn/sumOfWeights());
      scale(_hist_WplusMinus_MuAbseta, crossSection()/picobarn/sumOfWeights());
    }


    Histo1DPtr _hist_Wplus_MuAbseta;
    Histo1DPtr _hist_Wminus_MuAbseta;
    Histo1DPtr _hist_WplusMinus_MuAbseta;

  };


  RIVET_DECLARE_PLUGIN(CMS_2019_I1705068);

}
