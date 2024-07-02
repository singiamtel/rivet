// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> c l nu moments
  class BELLE_2023_I2649712 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2023_I2649712);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521),"UFS");
      // histos
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_p[ix  ], 1+ix, 1, 1);
        book(_p[ix+4], "TMP/p_"+toString(ix),refData<YODA::BinnedEstimate<string>>(1+ix, 1, 1));
      }
    }

    void findDecayProducts(Particle parent, Particles& em, Particles& ep,
                           Particles& nue, Particles& nueBar, bool& charm) {
      for (const Particle& p : parent.children()) {
        if (PID::isCharmHadron(p.pid())) {
          charm=true;
        }
        else if (p.pid() == PID::EMINUS || p.pid()==PID::MUON) {
          em.push_back(p);
        }
        else if (p.pid() == PID::EPLUS || p.pid()==PID::ANTIMUON) {
          ep.push_back(p);
        }
        else if (p.pid() == PID::NU_E  || p.pid()==PID::NU_MU) {
          nue.push_back(p);
        }
        else if (p.pid() == PID::NU_EBAR || p.pid()==PID::NU_MUBAR) {
          nueBar.push_back(p);
        }
        else if (PID::isBottomHadron(p.pid())) {
          findDecayProducts(p,em,ep,nue,nueBar,charm);
        }
        else if (!PID::isHadron(p.pid())) {
          findDecayProducts(p,em,ep,nue,nueBar,charm);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_edges.empty()) _edges=_p[0]->xEdges();
      // find and loop over Upslion(4S)
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        if (p.children().empty() || (p.children().size()==1 && p.children()[1].abspid()==p.abspid())) {
          continue;
        }
        // find decay products
        bool charm = false;
        Particles em,ep,nue,nueBar;
        findDecayProducts(p,em,ep,nue,nueBar,charm);
        if (!charm) continue;
        FourMomentum pl,pnu;
        if (em.size()==1 && nueBar.size()==1 && em[0].pid()+1==-nueBar[0].pid()) {
          pl  = em[0].mom();
          pnu = nueBar[0].mom();
        }
        else if(ep.size()==1 && nue.size()==1 && nue[0].pid()==-ep[0].pid()+1) {
          pl  = ep[0].mom();
          pnu = nue[0].mom();
       	}
        else {
          continue;
        }
        double q2 = (pl+pnu).mass2();
        vector<double> q2n(8);
        for (unsigned int ix=0;ix<8;++ix) q2n[ix] = pow(q2,1+ix);
        if (q2<1.5) continue;
        double q2cut=1.5;
        for (unsigned int ibin=0; ibin<15; ++ibin) {
          if (q2>q2cut) {
            for (unsigned int ix=0; ix<8; ++ix) {
              _p[ix]->fill(_edges[ibin],q2n[ix]);
            }
          }
          else break;
          q2cut+=0.5;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> tmp[3];
      for (unsigned int ix=0; ix<3;++ix) book(tmp[ix],5+ix,1,1);
      for (unsigned int iy=0; iy<_p[0]->numBins();++iy) {
        double q2  = _p[0]->bin(iy+1).mean(2)  ,  q4  = _p[1]->bin(iy+1).mean(2),
               q6  = _p[2]->bin(iy+1).mean(2)  ,  q8  = _p[3]->bin(iy+1).mean(2),
              q10  = _p[4]->bin(iy+1).mean(2)  ,  q12 = _p[5]->bin(iy+1).mean(2),
          q14  = _p[6]->bin(iy+1).mean(2)  ,  q16 = _p[7]->bin(iy+1).mean(2);
        // q4 case
        double N = _p[0]->bin(iy+1).effNumEntries();
        double value = q4 - sqr(q2);
        double error = (-sqr(q4) + 4*sqr(q2)*(-sqr(q2) + 2*q4) - 4*q2*q6 + q8)/N;
        tmp[0]->bin(iy+1).set(value,sqrt(error));
        // q6 case
        value = q6 + q2*(2*sqr(q2) - 3*q4);
        error = (q12 - sqr(q6) - 6*q2*(q10 - 5*q4*q6) + 3*q4*(3*sqr(q4) - 2*q8) + 
                 sqr(q2)*(-72*sqr(q4) + 36*sqr(q2)*(-sqr(q2) + 3*q4) - 48*q2*q6 + 21*q8))/N;
        tmp[1]->bin(iy+1).set(value,sqrt(error));
        // q8 case
        value = q8 + q2*(-3*q2*sqr(q2) + 6*q2*q4 - 4*q6);
        error = (q16 - 8*q6*(q10 - 2*q4*q6) - sqr(q8) - 8*q2*(q14 - 3*q4*(q10 - 4*q4*q6) - 6*q6*q8) + 
                 4*sqr(q2)*(7*q12 - 28*sqr(q6) + 6*q2*(-3*q10 + 22*q4*q6) + 3*q4*(12*sqr(q4) - 11*q8) + 
                            3*sqr(q2)*(-12*sqr(q2)*(sqr(q2) - 4*q4) - 51*sqr(q4) - 28*q2*q6 + 13*q8)))/N;
        tmp[2]->bin(iy+1).set(value,sqrt(error));
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedProfilePtr<string> _p[8];
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2023_I2649712);

}
