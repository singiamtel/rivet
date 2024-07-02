// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+e-> Z+- pi-+
  class BESIII_2015_I1391798 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1391798);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9044213);
      // projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      // histograms
      unsigned int iEnergy=0;
      if (isCompatibleWithSqrtS(4.23)) {
        _ecms="4.23";
        iEnergy=1;
      }
      else if (isCompatibleWithSqrtS(4.26)) {
        _ecms="4.26";
        iEnergy=2;
      }
      else {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
      // histos
      for (unsigned int ix=0;ix<2;++ix) {
        book(_h[ix],2,ix+1,iEnergy);
        book(_h[2+ix],3,1,1+ix);
        book(_sigma[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.pid()>0) {
	      axis = beams.first .momentum().p3().unit();
      }
      else {
        axis = beams.second.momentum().p3().unit();
      }
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DD,other;
      for (const Particle& p : fs) {
        Particle parent=p;
        while (!parent.parents().empty()) {
          parent=parent.parents()[0];
          if (parent.abspid()==411 || parent.abspid()==413 ||
	          parent.abspid()==421 || parent.abspid()==423) break;
        }
        if ((parent.abspid()==411 || parent.abspid()==421)
           && !parent.parents().empty()) {
          Particle Dstar = parent.parents()[0];
          if (Dstar.abspid()==413 || Dstar.abspid()==423) {
            parent=Dstar;
          }
        }
        if (parent.abspid()==411 || parent.abspid()==413 ||
           parent.abspid()==421 || parent.abspid()==423) {
          bool found=false;
          for (const auto& D : DD) {
            // D already in list
            if (fuzzyEquals(D.mom(), parent.mom())) {
              found=true;
              break;
            }
          }
          if(!found) DD.push_back(parent);
        }
        else {
          other.push_back(p);
        }
      }
      // D Dbar + charged pion
      if (DD.size()!=2 || other.size()!=1) vetoEvent;
      if (DD[0].pid()*DD[1].pid()>0) vetoEvent;
      if (other[0].abspid()!=211) vetoEvent;
      if (DD[0].abspid()%10!=3) swap(DD[0],DD[1]);
      if (DD[0].abspid()%10!=3 || DD[1].abspid()%10!=1) vetoEvent;
      const double mass = (DD[0].momentum()+DD[1].momentum()).mass();
      const double cTheta = abs(axis.dot(other[0].momentum().p3().unit()));
      unsigned int iloc=0;
      // D0 D*- pi+ +cc
      if ((DD[0].pid()==-413 && DD[1].pid()== 421 && other[0].pid()== 211) ||
      	 (DD[0].pid()== 413 && DD[1].pid()==-421 && other[0].pid()==-211)) {
       	iloc=0;
      }
      // D- D*0 pi+ +cc
      else if((DD[0].pid()== 423 && DD[1].pid()==-411 && other[0].pid()== 211) ||
	       (DD[0].pid()==-423 && DD[1].pid()== 411 && other[0].pid()==-211)) {
      	iloc=1;
      }
      // otherwise veto event
      else {
      	vetoEvent;
      }
      _h[iloc]->fill(mass);
      // parent Z+/-
      if (DD[0].parents()[0].abspid()==_pid && DD[1].parents()[0].abspid()==_pid &&
      	 fuzzyEquals(DD[0].parents()[0].momentum(),DD[1].parents()[0].momentum()) ) {
        _sigma[0]->fill(_ecms);
        if (_ecms=="4.26") _sigma[1]->fill(_ecms);
        _h[2+iloc]->fill(cTheta);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // distributions
      normalize(_h, 1.0, false);
      // cross section
      for(unsigned int ix=0;ix<2;++ix)
        scale(_sigma[ix], crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    Histo1DPtr _h[4];
    int _pid;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1391798);

}
