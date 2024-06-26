// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- -> double charmonium
  class BELLE_2004_I653673 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I653673);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare("FS",FinalState());
      declare("UFS",UnstableParticles(Cuts::pid==443 || Cuts::pid==100443 ||
                    Cuts::pid==441 || Cuts::pid==10441 || Cuts::pid==100441));
      // histograms
      for (unsigned int ix=0; ix<3; ++ix) {
        for (unsigned int iy=0;iy<3;++iy) {
          book(_p[ix][iy],"TMP/p_"+toString(ix+1)+"_"+toString(iy+1),
               refData<YODA::BinnedEstimate<string>>(3,1+ix,1+iy));
          if (ix==2) continue;
          book(_h_sigma[ix][iy], 1+ix, 1   , 1+2*iy);
          book(_h_angle[ix][iy], 4   , 1+ix, 1+iy  );
        }
      }
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int& ncount,
                      unsigned int & nCharged) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
          if (PID::isCharged(p.pid())) ++nCharged;
        }
        else {
          findChildren(child,nRes,ncount,nCharged);
        }
      }
    }

    double helicityAngle(const Particle & p) const {
      if (p.children().size()!=2) return 10.;
      if (p.children()[0].abspid()!=PID::MUON || p.children()[0].pid()!=-p.children()[1].pid()) {
        return 10.;
      }
      Vector3 axis = p.p3().unit();
      const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
      for (const Particle& child : p.children()) {
        if (child.pid()!=PID::MUON) continue;
        return axis.dot(boost.transform(child.momentum()).p3().unit());
      }
      return 10;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // loop over J/psi and psi(2S)
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      bool matched=false;
      for (const Particle& p : ufs.particles(Cuts::pid==443 || Cuts::pid==100443 )) {
      	if (p.children().empty()) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
        unsigned int nCharged=0;
      	findChildren(p,nRes,ncount,nCharged);
        // eta_c, chi_c0, eta_c(2S)
      	for (const Particle& p2 : ufs.particles(Cuts::pid==441 || Cuts::pid==10441 || Cuts::pid==100441)) {
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          unsigned int nChargedB=0;
          findChildren(p2,nResB,ncountB,nChargedB);
          if (ncountB!=0) continue;
          matched = true;
          for (const auto& val : nResB) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            unsigned int ipsi = p.pid()/100000;
            unsigned int ieta = p2.pid()/10000;
            if (ieta>1) ieta=2;
            // fill the cross sections
            if ((ipsi==0 && nChargedB>2) || (ipsi==1 && nChargedB>0)) {
              _h_sigma[ipsi][ieta]->fill(_ecms);
            }
            if (ipsi>0) break;
            // angular dists for J/psi only
            // production
            const double cProd = p.p3().z()/p.p3().mod();
            _h_angle[0][ieta]->fill(abs(cProd));
            _p[0][ieta]->fill(_ecms,-1.25*(1.-3.*sqr(cProd)));
            _p[2][ieta]->fill(_ecms,-1.25*(1.-3.*sqr(cProd)));
            // helicity angle
            const double cHel = helicityAngle(p);
            if (cHel>1.) break;
            _h_angle[1][ieta]->fill(abs(cHel));
            _p[1][ieta]->fill(_ecms,-1.25*(1.-3.*sqr(cHel)));
            _p[2][ieta]->fill(_ecms,-1.25*(1.-3.*sqr(cHel)));
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_sigma[ix], crossSection()/sumOfWeights()/femtobarn);
        normalize(_h_angle[ix]);
      }
      // extract the alpha parameters
      for (unsigned int ix=0; ix<3; ++ix) {
        for (unsigned int iy=0; iy<3; ++iy) {
          const double val = _p[ix][iy]->bin(1).mean(2);
          const double err = _p[ix][iy]->bin(1).stdErr(2);
          BinnedEstimatePtr<string> tmp;
          book(tmp,3,1+ix,1+iy);
          const double alpha = 3.*val/(1-val);
          const double error = 3./sqr(1.-val)*err;
          tmp->bin(1).set(alpha,error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_sigma[2][3];
    Histo1DPtr _h_angle[2][3];
    BinnedProfilePtr<string> _p[3][3];
    string _ecms="10.6";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2004_I653673);

}
