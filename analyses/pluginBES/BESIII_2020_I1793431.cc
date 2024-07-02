// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- -> pi0 pi0 J/psi
  class BESIII_2020_I1793431 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1793431);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      // histos
      unsigned int iloc=0;
      if      (isCompatibleWithSqrtS(4.226)) iloc=1;
      else if (isCompatibleWithSqrtS(4.236)) iloc=2;
      else if (isCompatibleWithSqrtS(4.244)) iloc=3;
      else if (isCompatibleWithSqrtS(4.258)) iloc=4;
      if (iloc>0) {
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_h[iy], 2, iloc, 1+iy);
        }
      }
      book(_sigma[0], 1, 1, 1);
      book(_sigma[1],"TMP/numer", refData<YODA::BinnedEstimate<string>>(1,1,2));
      book(_sigma[2],"TMP/denom", refData<YODA::BinnedEstimate<string>>(1,1,2));
      book(_sigma[3], 3, 1, 1);

      for (unsigned int ix=0; ix<3; ++ix) {
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles Jpsi, other;
      for (const Particle& p : fs) {
        Particle parent=p;
        while (!parent.parents().empty()) {
          parent=parent.parents()[0];
          if (parent.abspid()==PID::JPSI) break;
        }
        if (parent.abspid()!=PID::JPSI) {
          other.push_back(p);
          continue;
        }
        bool found=false;
        for (auto & psi : Jpsi) {
          // J/psi already in list
          if (fuzzyEquals(psi.momentum(),parent.momentum())) {
            found=true;
            break;
          }
        }
        if (!found) Jpsi.push_back(parent);
      }
      if (Jpsi.size()!=1 || other.size()!=2) vetoEvent;
      if (other[0].pid()==PID::PI0 && other[1].pid()==PID::PI0) {
        _sigma[0]->fill(_ecms[0]);
        if(!_ecms[1].empty()) _sigma[1]->fill(_ecms[1]);
        if (Jpsi[0].parents()[0].pid()==_pid && !_ecms[2].empty()) _sigma[3]->fill(_ecms[2]);
        if (_h[0]) {
          for (unsigned int iy=0; iy<2; ++iy) {
            _h[0]->fill((Jpsi[0].mom()+other[iy].mom()).mass());
          }
          _h[1]->fill((other[0].mom()+other[1].mom()).mass());
        }
      }
      else if (other[0].pid()==-other[1].pid() &&
               other[0].abspid()==PID::PIPLUS &&
               !_ecms[1].empty()) {
        _sigma[2]->fill(_ecms[1]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int iy=0; iy<2; ++iy) {
	      if (_h[iy]) normalize(_h[iy], 1.0, false);
      }
      scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
      BinnedEstimatePtr<string> tmp;
      book(tmp,1,1,2);
      divide(_sigma[1], _sigma[2], tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h[2];
    BinnedHistoPtr<string> _sigma[4];
    string _ecms[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1793431);

}
