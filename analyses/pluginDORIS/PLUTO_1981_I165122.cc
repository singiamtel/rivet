// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief kaon production at low energies
  class PLUTO_1981_I165122 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1981_I165122);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // // Book histograms
      book(_c_hadrons , "/TMP/sigma_hadrons", refData<YODA::BinnedEstimate<string>>(1, 1, 1));
      for(unsigned int ix=0;ix<2;++ix) {
        book(_c_muons[ix], "/TMP/sigma_muons_"+toString(ix), refData<YODA::BinnedEstimate<string>>(1+2*ix, 1, 1));
        book(_c_kaons[ix], "/TMP/sigma_kaons_"+toString(ix), refData<YODA::BinnedEstimate<string>>(1+2*ix, 1, 1));
        for (const string& en : _c_muons[ix].binning().edges<0>()) {
          // these points are really Upsilon(1S) decay
          if(en == "9.45 - 9.466"s || en=="9.45 - 9.456"s) continue;
          const size_t idx = en.find("-");
          if(idx!=std::string::npos) {
            const double emin = std::stod(en.substr(0,idx));
            const double emax = std::stod(en.substr(idx+1,string::npos));
            if(inRange(sqrtS()/GeV, emin, emax)) {
              _ecms[ix] = en;
              break;
            }
          }
          else {
            const double end = std::stod(en)*GeV;
            if (isCompatibleWithSqrtS(end)) {
              _ecms[ix] = en;
              break;
            }
          }
        }
      }
      if (isCompatibleWithSqrtS(9.4*GeV)) {
        book(_h_spectrum1, 5, 1, 1);
      }
      else if (isCompatibleWithSqrtS(30.0*GeV, 1E-2)) {
        book(_h_spectrum1, 4, 1, 1);
      }
      book(_h_spectrum2, 6, 1, 1);
      book(_c_hadronsY,"TMP/nUps");
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles &kaons, Particles& stable) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
        if(id==130 || id ==310) {
          kaons.push_back(p);
        }
        if (id==111 or p.children().empty())
          stable.push_back(p);
        else
          findDecayProducts(p, kaons, stable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      // Continuum
      if (upsilons.empty()) {
        MSG_DEBUG("No Upsilons found => continuum event");
        // final state particles
        const FinalState& fs = apply<FinalState>(event, "FS");
        map<long,int> nCount;
        int ntotal(0);
        for (const Particle& p : fs.particles()) {
          nCount[p.pid()] += 1;
          ++ntotal;
        }
        // mu+mu- + photons
        if(nCount[-13]==1 and nCount[13]==1 &&
           ntotal==2+nCount[22]) {
          for(unsigned int ix=0;ix<2;++ix) _c_muons[ix]->fill(_ecms[ix]);
        }
        // everything else
        else
          _c_hadrons->fill(_ecms[0]);
        // unstable particles
        for (const Particle& p : ufs.particles(Cuts::pid==130 or Cuts::pid==310)) {
          if(_h_spectrum1) {
            double xp = p.p3().mod()/meanBeamMom;
            _h_spectrum1->fill(xp);
          }
          for(unsigned int ix=0;ix<2;++ix) _c_kaons[ix]->fill(_ecms[ix]);
        }
      }
      else {
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
          Particles kaons,stable;
          // Find the decay products we want
          findDecayProducts(ups, kaons, stable);
          // boost to rest frame (if required)
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();

          map<long,int> nCount;
          int ntotal(0);
          for (const Particle& p : stable) {
            nCount[p.pid()] += 1;
            ++ntotal;
          }
          for( const Particle & kaon : kaons) {
            const FourMomentum p2 = cms_boost.transform(kaon.momentum());
            const double xp = 2.*p2.p3().mod()/mass;
            _h_spectrum2->fill(xp);
            _c_kaons[0]->fill("9.45 - 9.466"s);
            _c_kaons[1]->fill("9.45 - 9.456"s);
          }
          // mu+mu- + photons
          if(nCount[-13]==1 and nCount[13]==1 &&
             ntotal==2+nCount[22]) {
            _c_muons[0]->fill("9.45 - 9.466"s);
            _c_muons[1]->fill("9.45 - 9.456"s);
          }
          // everything else
          else {
            _c_hadrons->fill("9.45 - 9.466"s);
            _c_hadronsY->fill();
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> ratio;
      book(ratio,1,1,1);
      divide(_c_kaons[0],_c_muons[0],ratio);
      book(ratio,2,1,1);
      divide(_c_kaons[0],_c_hadrons,ratio);
      book(ratio,3,1,1);
      divide(_c_kaons[1],_c_muons[1],ratio);
      // normalize the spectra if required
      if(_h_spectrum1) {
        scale(_h_spectrum1, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      }
      if(_h_spectrum2) {
        scale(_h_spectrum2, 1./ *_c_hadronsY);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons,_c_muons[2],_c_kaons[2];
    string _ecms[2];
    Histo1DPtr _h_spectrum1, _h_spectrum2;
    CounterPtr _c_hadronsY;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1981_I165122);


}
