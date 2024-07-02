// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xc ell - nu_ell
  class BELLE_2007_I732595 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I732595);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521),"UFS");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_p[ix], 1, 1, 1+2*ix);
      }
    }

    void findDecayProducts(const Particle& parent, Particles& em, Particles& ep,
                           Particles& nue, Particles& nueBar, bool& charm) {
      for (const Particle & p : parent.children()) {
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
      if(_edges.empty()) {
        _edges = _p[0]->xEdges();
        for(const string & en : _edges)
          _eCut.push_back(std::stod(en)*GeV);
      }
      // find and loop over B mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        if (p.children().empty() ||
           (p.children().size()==1 && p.children()[1].abspid()==p.abspid())) continue;
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
        else if (ep.size()==1 && nue.size()==1 && nue[0].pid()==-ep[0].pid()+1) {
          pl  = ep[0].mom();
          pnu = nue[0].mom();
        }
        else {
          continue;
        }
        // boost to rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        double mX2 = (p.mom()-pl-pnu).mass2();
        pl = boost.transform(pl);
        double Estar = pl.E();
        for (unsigned int ix=0;ix<2;++ix) {
          for(unsigned int ix=0;ix<_eCut.size();++ix) {
            if(Estar>_eCut[ix]) {
              _p[0]->fill(_edges[ix],mX2);
              _p[1]->fill(_edges[ix],sqr(mX2));
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // compute <(mx2-<mx2>)^2> = <mx4>-<mx2>^2
      BinnedEstimatePtr<string> tmp;
      book(tmp,1,1,2);
      for (unsigned int ix=0; ix<_p[0]->numBins(); ++ix) {
        const double value = _p[1]->bin(ix+1).mean(2)-sqr(_p[0]->bin(ix+1).mean(2));
        const double error = value*sqrt(sqr(_p[1]->bin(ix+1).relStdErr(2))+4.*sqr(_p[0]->bin(ix+1).relStdErr(2)));
        tmp->bin(ix+1).set(value,error);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedProfilePtr<string> _p[2];
    vector<string> _edges;
    vector<double> _eCut;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I732595);

}
