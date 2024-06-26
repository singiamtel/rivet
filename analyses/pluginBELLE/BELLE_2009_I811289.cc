// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief BELLE J/psi production
  class BELLE_2009_I811289 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I811289);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare("FS",FinalState());
      declare("CFS",FinalState());
      // book histos
      for (unsigned int ix=0; ix<5; ++ix) {
        book(_h_spect[ix], 2, 1, 1+ix);
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_sigma[ix], 1, 1, 1+ix);
        for (unsigned int iy=0; iy<3; ++iy) {
          book(_h_angle[ix][iy], 3, 1+ix, 1+iy);
        }
      }
    }

    void plotHelicityAngle(const Particle& p, Histo1DPtr h) {
      if (p.children().size()!=2) return;
      if (p.children()[0].abspid()!=PID::MUON || p.children()[0].pid()!=-p.children()[1].pid()) return;
      Vector3 axis = p.p3().unit();
      const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
      for (const Particle& child : p.children()) {
        if (child.pid()!=PID::MUON) continue;
        h->fill(axis.dot(boost.transform(child.mom()).p3().unit()));
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // 4 charged particle veto
      if (apply<FinalState>(event,"CFS").particles().size()<4) vetoEvent;
      // find the charmonium
      Particles D,Jpsi_direct,Jpsi_feed,onium;
      for (const Particle & p : apply<FinalState>(event,"FS").particles()) {
        Particle parent = p.parents()[0], Jpsi;
        while (!parent.parents().empty()) {
          // charm mesons
          if (parent.abspid()==411 || parent.abspid()==421 || parent.abspid()==431) {
            break;
          }
          // c cbar mesons
          if ((parent.abspid()%1000)/10 == 44 ) {
            if ((parent.parents()[0].abspid()%1000)/10!=44) break;
            if (parent.pid()==443) Jpsi=parent;
          }
          parent = parent.parents()[0];
        }
        // charm
        if (parent.abspid()==411 || parent.abspid()==421 || parent.abspid()==431) {
          bool found=false;
          for (auto & p : D) {
            // D already in list
            if (fuzzyEquals(p.mom(),parent.mom())) {
              found=true;
              break;
            }
          }
          if (!found) D.push_back(parent);
        }
        // c cbar
        else if ((parent.abspid()%1000)/10 == 44) {
          if (parent.pid()==443) {
            bool found=false;
            for (const auto& p : Jpsi_direct) {
              // psi already in list
              if (fuzzyEquals(p.mom(),parent.mom())) {
                found=true;
                break;
              }
            }
            if (!found) Jpsi_direct.push_back(parent);
          }
          else {
            bool found=false;
            for (const auto & p : onium) {
              // psi already in list
              if (fuzzyEquals(p.mom(),parent.mom())) {
                found=true;
                break;
              }
            }
            if (!found) onium.push_back(parent);
          }
        }
        if (Jpsi.pid()==443) {
          bool found=false;
          for (const auto& p : Jpsi_feed) {
            // psi already in list
            if (fuzzyEquals(p.mom(),Jpsi.mom())) {
              found=true;
              break;
            }
          }
          if (!found) Jpsi_feed.push_back(Jpsi);
        }
      }
      if (Jpsi_feed.empty() && Jpsi_direct.empty()) vetoEvent;
      // cross sections first
      if (D.empty() && onium.empty() && Jpsi_direct.size()==1) {
        _h_sigma[1]->fill("10.6"s);
        for (const Particle & p : Jpsi_direct) {
          double modp=p.p3().mod();
          _h_spect[4]->fill(modp);
          plotHelicityAngle(p,_h_angle[0][2]);
          _h_angle[1][2]->fill(abs(p.p3().z()/modp));
        }
      }
      else {
        _h_sigma[0]->fill("10.6"s);
        for (const Particle& p : Jpsi_direct) {
          double modp=p.p3().mod();
          _h_spect[3]->fill(modp);
          plotHelicityAngle(p,_h_angle[0][1]);
          _h_angle[1][1]->fill(abs(p.p3().z()/modp));
        }
        for (const Particle& p : Jpsi_feed) {
          double modp=p.p3().mod();
          _h_spect[3]->fill(modp);
          plotHelicityAngle(p,_h_angle[0][1]);
          _h_angle[1][1]->fill(abs(p.p3().z()/modp));
        }
        if (!D.empty()) {
          for (const Particle& p : Jpsi_direct) {
            _h_spect[1]->fill(p.p3().mod());
          }
          for (const Particle& p : Jpsi_feed) {
            _h_spect[1]->fill(p.p3().mod());
          }
        }
        else {
          for (const Particle& p : Jpsi_direct) {
            _h_spect[2]->fill(p.p3().mod());
          }
          for (const Particle& p : Jpsi_feed) {
            _h_spect[2]->fill(p.p3().mod());
          }
        }
      }
      for (const Particle& p : Jpsi_direct) {
        double modp=p.p3().mod();
        _h_spect[0]->fill(modp);
        plotHelicityAngle(p,_h_angle[0][0]);
        _h_angle[1][0]->fill(abs(p.p3().z()/modp));
      }
      for (const Particle& p : Jpsi_feed) {
        double modp=p.p3().mod();
        _h_spect[0]->fill(modp);
        plotHelicityAngle(p,_h_angle[0][0]);
        _h_angle[1][0]->fill(abs(p.p3().z()/modp));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spect, crossSection()/sumOfWeights()/femtobarn);
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_sigma[ix], crossSection()/sumOfWeights()/picobarn);
        normalize(_h_angle[ix], 1.0, false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_sigma[2];
    Histo1DPtr _h_spect[5],_h_angle[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I811289);

}
