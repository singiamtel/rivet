// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {


  /// @brief D spect in B decays
  class CLEOII_1997_I440969 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1997_I440969);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553),"UFS");
      // histos
      const vector<double> xbins={0.0,0.1,0.2,0.3,0.4,0.5};
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
      book(_b[0], {0.1,0.2,0.3,0.4,0.5});
      book(_b[1], {0.0,0.1,0.2,0.3,0.4,0.5});
      for (auto& b : _b[1]->bins()) {
        book(b, 6, 1, b.index());
        if (b.index() < 5) {
          book(_b[0]->bin(b.index()), 5, 1, b.index());
        }
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_p[ix],"TMP/p_"+toString(ix+1), refData(7,1,1+ix));
      }
    }

    void findDecay(const Particle& parent, Particles & Ds) {
      for (const Particle& p : parent.children()) {
        if (p.abspid()==421 || p.abspid()==411) {
          Ds.push_back(p);
        }
        else if (p.abspid()==413 || p.abspid()==423) {
          Ds.push_back(p);
          findDecay(p,Ds);
        }
        else if (p.abspid()==PID::PIPLUS ||
          p.abspid()==PID::KPLUS  ||
          p.pid()==PID::PI0 ||
          p.pid()==PID::K0S) {
          continue;
        }
        else if (!p.children().empty()) {
          findDecay(p,Ds);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        for (const Particle& B:p.children()) {
          if (B.abspid()!=511 && B.abspid()!=521) continue;
          Particles Ds;
          findDecay(B,Ds);
          for (const Particle & p2 : Ds) {
            FourMomentum pD = boost.transform(p2.mom());
            const double x = pD.p3().mod()/sqrt(0.25*sqr(p.mass())-sqr(p2.mass()));
            if (p2.abspid()==421)      _h[0]->fill(x);
            else if (p2.abspid()==411) _h[1]->fill(x);
            else {
              // x distribution for D*
              if (p2.abspid()==423) _h[2]->fill(x);
              else                  _h[3]->fill(x);
              // find decay products
              if (p2.children().size()!=2) continue;
              Particle dec;
              if ((p2.children()[0].abspid()==411 || p2.children()[0].abspid()==421) &&
                  (p2.children()[1].abspid()==111 || p2.children()[1].abspid()==211)) {
                dec = p2.children()[0];
              }
              else if((p2.children()[1].abspid()==411 || p2.children()[1].abspid()==421) &&
                      (p2.children()[0].abspid()==111 || p2.children()[0].abspid()==211)) {
                dec = p2.children()[1];
              }
              else continue;
              const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pD.betaVec());
              FourMomentum pDec = boost2.transform(boost.transform(dec.momentum()));
              // decay angle
              double cTheta = pDec.p3().unit().dot(pD.p3().unit());
              if (p2.abspid()==423) {
                _b[1]->fill(x,cTheta);
                _p[1]->fill(x,-1./2.*(1.-5.*sqr(cTheta)));
              }
              else {
                _b[0]->fill(x,cTheta);
                _p[0]->fill(x,-1./2.*(1.-5.*sqr(cTheta)));
              }
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // spectra
      normalize(_h, 1.0, false);
      for (unsigned int ix=0; ix<2; ++ix) {
        // angles dists
        normalize(_b[ix], 1.0, false);
        // polarizations
        Estimate1DPtr tmp;
        book(tmp, 7, 1, 1+ix);
        for (const auto& b : _p[ix]->bins()) {
          const double val = b.numEntries()>0 && b.effNumEntries()>0? b.yMean() : 0.;
          const double err = b.numEntries()>1 && b.effNumEntries()>1? b.yStdErr() : 0.;
          const double l1 = 2.*val/(1.+val);
          const double e1 = 2./sqr(1.+val)*err;
          tmp->bin(b.index()).set(l1, e1);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    Profile1DPtr _p[2];
    Histo1DGroupPtr _b[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1997_I440969);

}
