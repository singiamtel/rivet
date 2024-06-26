// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> c l nu moments
  class CLEOII_2004_I647287 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2004_I647287);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521),"UFS");
      // spin weighted average D and D* mass
      _mD2 = sqr(0.125*(1.86966+1.86484)+0.375*(2.00685+2.01026));
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<4; ++iy) {
          if (iy%2==1) {
            book(_p[ix][iy],
                 "TMP/p_"+toString(ix)+"_"+toString(iy),
                 refData<YODA::BinnedEstimate<string>>(1,1+iy,1+ix));
          }
          else {
            book(_p[ix][iy],1,+1+iy,1+ix);
          }
        }
      }
      book(_p_dist,2,1,1);
    }

    void findDecayProducts(const Particle& parent, Particles& em, Particles& ep,
                           Particles& nue, Particles& nueBar, bool & charm) {
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
      if(_edges.empty()) _edges=_p_dist->xEdges();
      const double Ecut[2]={1.,1.5};
      // find and loop over Upslion(4S)
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        if (p.children().empty() ||
           (p.children().size()==1 && p.children()[1].abspid()==p.abspid())) {
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
        else if (ep.size()==1 && nue.size()==1 && nue[0].pid()==-ep[0].pid()+1) {
          pl  = ep[0].mom();
          pnu = nue[0].mom();
       	}
        else {
          continue;
        }
        // boost to rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        double q2 = (pl+pnu).mass2();
        FourMomentum pX = boost.transform(p.mom()-pl-pnu);
        pl = boost.transform(pl);
        double mX2  = pX.mass2();
        for (unsigned int ix=0;ix<2;++ix) {
          if (pl.E()>Ecut[ix]) {
            _p[ix][0]->fill(_edges[5*ix],mX2-_mD2);
            _p[ix][1]->fill(_edges[5*ix],sqr(mX2));
            _p[ix][2]->fill(_edges[5*ix],q2);
            _p[ix][3]->fill(_edges[5*ix],sqr(q2));
          }
        }
        for(unsigned int ix=0;ix<6;++ix) {
          if (1.+0.1*ix < pl.E())  _p_dist->fill(_edges[ix],mX2-_mD2);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        // compute <(mx2-<mx2>)^2> = <mx4>-<mx2>^2
        BinnedEstimatePtr<string> tmp;
        book(tmp,1,2,1+ix);
        for (const auto& b0 : _p[ix][0]->bins()) {
          const auto& b1 = _p[ix][1]->bin(b0.index());
          const double val = b1.yMean()-sqr(b0.yMean()+_mD2);
          const double err = val*sqrt(sqr(b1.relErrW())+4.*sqr(b0.yStdErr()/(b0.yMean()+_mD2)));
          tmp->bin(b0.index()).set(val, err);
        }
        // compute <(q2-<q2>)^2> = <q4>-<q2>^2
        book(tmp,1,4,1+ix);
        for (const auto& b2 : _p[ix][2]->bins()) {
          const auto& b3 = _p[ix][3]->bin(b2.index());
          const double val = b3.yMean()-sqr(b2.yMean());
          const double err = val*sqrt(sqr(b3.relErrW())+4.*sqr(b3.relErrW()));
          tmp->bin(b2.index()).set(val, err);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedProfilePtr<string> _p[2][4], _p_dist;
    vector<string> _edges;
    double _mD2;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2004_I647287);

}
