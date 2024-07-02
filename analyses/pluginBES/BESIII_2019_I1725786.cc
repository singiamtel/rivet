// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > D D pi pi
  class BESIII_2019_I1725786 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1725786);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      // histograms
      unsigned int iloc=0;
      if      (isCompatibleWithSqrtS(4.4156)) {
        iloc=2;
        for (unsigned int ix=0; ix<2; ++ix) {
          book(_h_Dpi[ix],1,1,1+ix);
        }
        book (_h_DD,2,1,3);
      }
      else if (isCompatibleWithSqrtS(4.258)) {
        book(_h_DD,2,1,1);
      }
      else if (isCompatibleWithSqrtS(4.3583)) {
        book(_h_DD,2,1,2);
        iloc=1;
      }
      else if (isCompatibleWithSqrtS(4.5995)) {
        iloc=3;
      }
      if (iloc>0) {
        for (unsigned int iy=0;iy<3;++iy) {
          book(_h_mass[iy],3,iloc,1+iy);
        }
      }
      iloc=0;
      for (unsigned int ix=0;ix<5;++ix) {
        if (ix==1)  continue;
        book(_sigma[ix], 4+ix,1,1);
        if (ix!=0 && ix!=2 )  continue;
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[iloc] = en;
            break;
          }
        }
        iloc+=1;
      }
      if(_ecms[0].empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DD,other;
      for (const Particle & p : fs) {
        Particle parent=p;
        while(!parent.parents().empty()) {
          parent=parent.parents()[0];
          if(parent.abspid()==PID::D0 ||
             parent.abspid()==PID::DPLUS) break;
        }
        if(parent.abspid()!=PID::D0 &&
           parent.abspid()!=PID::DPLUS) {
          other.push_back(p);
          continue;
        }
        bool found=false;
        for (auto& D : DD) {
          // D already in list
          if (fuzzyEquals(D.momentum(),parent.momentum())) {
            found=true;
            break;
          }
        }
        if(!found) DD.push_back(parent);
      }
      //  d dbar + 2 other particles
      if(DD.size()!=2 || other.size()!=2) vetoEvent;
      // other particles pi+ pi-
      if(!(other[0].pid()==-other[1].pid() && other[0].abspid()==PID::PIPLUS)) vetoEvent;
      // d dbar pair
      if(DD[0].pid()!=-DD[1].pid()) vetoEvent;
      if(other[0].pid()<0) swap(other[0],other[1]);
      if(DD   [0].pid()<0) swap(DD   [0],DD   [1]);
      // fill the mass plots if needed
      if (_h_Dpi[0] && DD[0].pid()==421) {
        _h_Dpi[0]->fill((DD[0].momentum()+other[0].momentum()).mass());
        _h_Dpi[1]->fill((DD[1].momentum()+other[1].momentum()).mass());
      }
      if (_h_DD) {
        _h_DD->fill((DD[0].momentum()+DD[1].momentum()).mass());
      }
      const bool DstarP = (DD[0].pid()== 421 && DD[0].parents()[0].pid()== 413);
      const bool DstarM = (DD[1].pid()==-421 && DD[1].parents()[0].pid()==-413);
      const bool Dstar = DstarP || DstarM;
      // D0 pi+ pi- mass
      if (_h_mass[0] && DD[0].pid()==421&&!Dstar) {
        _h_mass[0]->fill((DD[0].momentum()+other[0].momentum()+other[1].momentum()).mass());
      }
      // D*+ pi- mass
      if (_h_mass[1] && DstarP) {
        _h_mass[1]->fill((DD[0].momentum()+other[0].momentum()+other[1].momentum()).mass());
      }
      // D+ pi+ pi- mass
      if (_h_mass[2] && DD[0].pid()==411) {
        _h_mass[2]->fill((DD[0].momentum()+other[0].momentum()+other[1].momentum()).mass());
      }
      // now for the counters for the cross sections
      if (DD[0].parents()[0].pid()==30443 && DD[1].parents()[0].pid()==30443) {
        _sigma[0]->fill(_ecms[0]);
      }
      else if (DD[0].pid()==421 && (DD[0].parents()[0].pid()==10423 || DD[1].parents()[0].pid()==-10423)) {
        _sigma[2]->fill(_ecms[1]);
      }
      else if (DD[0].pid()==421 && ((DstarP && DD[0].parents()[0].parents()[0].pid()== 10423) ||
				   (DstarM && DD[1].parents()[0].parents()[0].pid()==-10423))) {
        _sigma[3]->fill(_ecms[1]);
      }
      else if (DD[0].pid()==411 && (DD[0].parents()[0].pid()==10413 || DD[1].parents()[0].pid()==-10413)) {
        _sigma[4]->fill(_ecms[1]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // dists
      if (_h_Dpi[0])  normalize(_h_Dpi,      1.0, false);
      if (_h_DD)      normalize(_h_DD,       1.0, false);
      if (_h_mass[0]) {
        for(unsigned int ix=0;ix<3;++ix)
          normalize(_h_mass[ix], 1.0, false);
      }
      // cross sections
      const vector<double> brs = { 0.93, 1., 1., 0.677, 1.};
      const double fact = crossSection()/ sumOfWeights() /picobarn;;
      for (unsigned int ii=0; ii<5; ++ii) {
	      if (ii==1) continue;
        scale(_sigma[ii], fact/brs[ii]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Dpi[2],_h_DD,_h_mass[3];
    BinnedHistoPtr<string> _sigma[5];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1725786);

}

