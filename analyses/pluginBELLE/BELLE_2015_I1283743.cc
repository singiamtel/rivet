// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Upsilon(5S) -> pi+pi- Upsilon(1,2,3S)
  class BELLE_2015_I1283743 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1283743);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(Beam(), "Beams");
      // histograms
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_sigma[ix],ix+1,1,2);
      }
      for (unsigned int ix=0;ix<12;++ix) {
        // book(_h_angle_ups2[ix],14,1,1+ix);
        if (ix>=8) continue;
        book(_h_mass_ups2 [ix],12,1,1+ix);
        // book(_h_angle_ups3[ix],15,1,1+ix);
        if (ix>=6) continue;
        book(_h_mass_all [ix],10,1,1+ix);
        book(_h_mass_ups1[ix],11,1,1+ix);
        book(_h_mass_ups3[ix],13,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis1;
      if (beams.first.pid()>0) {
        axis1 = beams.second.mom().p3().unit();
      }
      else {
        axis1 = beams.first .mom().p3().unit();
      }
      // find upsilon pi+ pi-
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles UPS, other;
      for(const Particle & p : fs) {
       	Particle parent=p;
        Particle ups=p;
       	while (!parent.parents().empty()) {
      	  parent=parent.parents()[0];
      	  if (parent.pid()==553 ||
              parent.pid()==100553||
              parent.pid()==200553) ups=parent;
      	}
       	if (ups.pid()%1000!=553) {
      	  other += p;
      	  continue;
      	}
      	bool found=false;
      	for (const auto& p : UPS) {
      	  // psi already in list
      	  if (fuzzyEquals(p.mom(), ups.mom())) {
      	    found=true;
      	    break;
      	  }
      	}
      	if (!found) UPS += ups;
      }
      // upsilon + 2 other particles
      if (UPS.size()!=1 || other.size()!=2) vetoEvent;
      // other particles pi+ pi-
      if (!(other[0].pid()==-other[1].pid() && other[0].abspid()==PID::PIPLUS)) vetoEvent;
      // cross secrtion
      int iups = UPS[0].pid()/100000;
      _h_sigma[iups]->fill(string("10.865"));
      // fill the mass plots
      double mpipi = (other[0].mom()+other[1].mom()).mass();
      double mUpspi[2] = {(other[0].mom()+UPS[0].mom()).mass(),
                          (other[1].mom()+UPS[0].mom()).mass()};
      if (mUpspi[0]<mUpspi[1]) {
        swap(other[0],other[1]);
        swap(mUpspi[0],mUpspi[1]);
      }
      double mpipi2=sqr(mpipi);
      double mMax2=sqr(mUpspi[0]);
      _h_mass_all[iups+3]->fill(mpipi);
      if (iups==0) {
        if (mpipi2>0.2)                  _h_mass_all [0]->fill(mUpspi[0]);
        if (mpipi2>0.2&&mpipi2<1)        _h_mass_ups1[0]->fill(mUpspi[0]);
        else if (mpipi2>1 && mpipi2<1.5) _h_mass_ups1[1]->fill(mUpspi[0]);
        else if (mpipi2>1.5)             _h_mass_ups1[2]->fill(mUpspi[0]);
        if (mMax2<105)                   _h_mass_ups1[3]->fill(mpipi);
        else if (mMax2>105&&mMax2<110)   _h_mass_ups1[4]->fill(mpipi);
        else if (mMax2>110)              _h_mass_ups1[5]->fill(mpipi);
      }
      else if (iups==1) {
        if (mpipi2>0.14)                 _h_mass_all [1]->fill(mUpspi[0]);
        if (mpipi2>0.14&&mpipi2<0.3)     _h_mass_ups2[0]->fill(mUpspi[0]);
        else if (mpipi2>0.3&&mpipi2<0.5) _h_mass_ups2[1]->fill(mUpspi[0]);
        else if (mpipi2>0.5&&mpipi2<0.6) _h_mass_ups2[2]->fill(mUpspi[0]);
        else if (mpipi2>0.6)             _h_mass_ups2[3]->fill(mUpspi[0]);
        if (mMax2<110)                   _h_mass_ups2[4]->fill(mpipi);
        else if (mMax2>110&&mMax2<112)   _h_mass_ups2[5]->fill(mpipi);
        else if (mMax2>112&&mMax2<114)   _h_mass_ups2[6]->fill(mpipi);
        else if (mMax2>114)              _h_mass_ups2[7]->fill(mpipi);
      }
      else if (iups==2) {
        if (mpipi2>0.1)                   _h_mass_all [2]->fill(mUpspi[0]);
        if (mpipi2>0.1&&mpipi2<0.15)      _h_mass_ups3[0]->fill(mUpspi[0]);
        else if (mpipi2>0.15&& mpipi2<.2) _h_mass_ups3[1]->fill(mUpspi[0]);
        else if (mpipi2>0.2)              _h_mass_ups3[2]->fill(mUpspi[0]);
        if (mMax2<113)                    _h_mass_ups3[3]->fill(mpipi);
        else if (mMax2>113&&mMax2<114)    _h_mass_ups3[4]->fill(mpipi);
        else if (mMax2>114)               _h_mass_ups3[5]->fill(mpipi);
      }
      // // only angular stuff for 2s and 3s
      // if(iups==0) return;
      // double cPi = axis1.dot(other[1].mom().p3().unit());
      // int iRegion=-1;
      // if (mUpspi[0]>10.605 && mUpspi[1]<10.635)
      // 	iRegion=0;
      // else if(mUpspi[0]>10.645 && mUpspi[1]<10.675)
      // 	iRegion=1;
      // else if(mUpspi[0]<10.57) {
      // 	iRegion=2;
      // 	if(iups==2) return;
      // }
      // else
      // 	return;
      // if(iups==1) {
      // 	_h_angle_ups2[4*iRegion]->fill(cPi);
      // }
      // else if(iups==2) {
      // 	_h_angle_ups3[4*iRegion]->fill(cPi);
      // }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma, crossSection()/ sumOfWeights() /picobarn);
      //normalize(_h_angle_ups2, 1.0, false);
      normalize(_h_mass_ups2, 1.0, false);
      //normalize(_h_angle_ups3, 1.0, false);
      normalize(_h_mass_all, 1.0, false);
      normalize(_h_mass_ups1, 1.0, false);
      normalize(_h_mass_ups3, 1.0, false);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_sigma[3];
    Histo1DPtr _h_mass_all[6],_h_mass_ups1[6],
      _h_mass_ups2[8],_h_mass_ups3[6];
    // Histo1DPtr _h_angle_ups2[12],_h_angle_ups3[8];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1283743);

}
