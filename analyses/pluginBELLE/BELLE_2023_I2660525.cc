// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s D0 production in e+e-
  class BELLE_2023_I2660525 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2023_I2660525);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      _mode=0;
      if      ( getOption("MODE") == "DECAY" ) _mode=0;
      else if ( getOption("MODE") == "SIGMA" ) _mode=1;
      // Upsilon decays
      if (_mode==0) {
        for (unsigned int ix=0;ix<2;++ix) {
          book(_c[ix],"TMP/n_Ups_"+toString(ix+4));
          book(_h_x [0][ix],3,1,ix+3);
          book(_h_x [1][ix],3,1,ix+1);
          for (unsigned int iy=0;iy<2;++iy) {
            book(_h_br[ix][iy],1,1,2*ix+iy+1);
          }
        }
        book(_h_fs,1,1,5);
      }
      else {
        for (unsigned int ix=0;ix<4;++ix) {
          book(_sigma[ix],2,1,1+ix);
        }
        for (const string& en : _sigma[0].binning().edges<0>()) {
          double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms = en;
            break;
          }
        }
        if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
      }
    }

    void findDecayProducts(const Particle& p, bool& hasBs, Particles& d0, Particles& ds) {
      for (const Particle& child : p.children()) {
        hasBs |= child.abspid()==531;
        if (child.abspid()==431) {
          ds.push_back(child);
        }
        else if (child.abspid()==421) {
          d0.push_back(child);
        }
        else if (!child.children().empty()) {
          findDecayProducts(child,hasBs,d0,ds);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Upsilon(4,5S) decay
      if (_mode==0) {
        const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
        Particles upsilons = ufs.particles(Cuts::pid==300553 || Cuts::pid==400553 || Cuts::pid==9000553);
        for (const Particle& ups : upsilons) {
          unsigned int iups = ups.pid()==300553 ? 0 : 1;
          _c[iups]->fill();
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
          Particles ds,d0;
          bool hasBs=false;
          findDecayProducts(ups, hasBs, d0, ds);
          if(iups==1&&hasBs) _h_fs->fill();
          for (const Particle & p : ds) {
            FourMomentum p2 = cms_boost.transform(p.mom());
            double x = p2.p3().mod()/sqrt(0.25*ups.mass2()-p2.mass2());
            _h_br[iups][0]->fill();
            _h_x [iups][0]->fill(x  );

          }
          for (const Particle & p : d0) {
            FourMomentum p2 = cms_boost.transform(p.mom());
            double x = p2.p3().mod()/sqrt(0.25*ups.mass2()-p2.mass2());
            _h_br[iups][1]->fill();
            _h_x [iups][1]->fill(x);
          }
        }
      }
      else {
        const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
        Particles B  = ufs.particles(Cuts::abspid==511 || Cuts::abspid==521);
        Particles Bs = ufs.particles(Cuts::abspid==531);
        // require B mesons
        if(B.empty() && Bs.empty()) vetoEvent;
        if(!B.empty()) _sigma[3]->fill(_ecms);
        Particles Ds = ufs.particles(Cuts::abspid==431);
        for (unsigned int ix=0;ix<Ds.size();++ix) {
          _sigma[0]->fill(_ecms);
        }
        Particles D0 = ufs.particles(Cuts::abspid==421);
        for (unsigned int ix=0;ix<D0.size();++ix) {
          _sigma[1]->fill(_ecms);
        }
        if (!Bs.empty() && !Ds.empty()) {
          for (const Particle& p : Bs) {
            if (p.children().size()==1) continue;
            Particles ds,d0;
            bool hasBs=false;
            findDecayProducts(p, hasBs, d0, ds);
            for (unsigned int ix=0;ix<ds.size();++ix) _sigma[2]->fill(_ecms);
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      if (_mode==0) {
      	for (unsigned int ix=0;ix<2;++ix) {
      	  if (_c[ix]->numEntries()==0) continue;
          scale(_h_x[ix], 1.0/ *_c[ix]);
          // brs in % and at 4S /2 as br is for B decays
          if (ix==0) scale(_h_br[ix], 50  / *_c[ix]);
          else       scale(_h_br[ix], 100./ *_c[ix]);
      	  if(ix==1) scale(_h_fs,100./ *_c[1]);
      	}
      }
      else {
        scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    unsigned int _mode;
    Histo1DPtr _h_x[2][2];
    CounterPtr _h_br[2][2],_h_fs,_c[2];
    BinnedHistoPtr<string> _sigma[4];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2023_I2660525);

}
