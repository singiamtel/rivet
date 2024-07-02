// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of R 3.65 and 3.872 GeV
  class BESII_2006_I735496 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2006_I735496);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==100443 || Cuts::abspid==30443),"UFS");
      declare(FinalState(), "FS");
      // Book histograms
      for (unsigned int ix=0;ix<2;++ix) {
        book(_sigma_hadrons[ix], "/TMP/sigma_hadrons_"+toString(ix+1),refData<YODA::BinnedEstimate<string>>(1,1,1+ix));
      }
      book(_sigma_muons,   "/TMP/sigma_muons",refData<YODA::BinnedEstimate<string>>(1,1,1));

      for (const string& en : _sigma_muons.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // mu+mu- + photons
      if (nCount[-13]==1 and nCount[13]==1 &&
         ntotal==2+nCount[22])
        _sigma_muons->fill(_ecms);
      // everything else
      else {
        _sigma_hadrons[1]->fill(_ecms);
        Particles psi = apply<UnstableParticles>(event,"UFS").particles();
        if (psi.empty()) {
          _sigma_hadrons[0]->fill(_ecms);
        }
        else {
          bool psi3770 = false;
          for (const Particle& p : psi) psi3770 |= p.pid()==30443;
          if (psi3770) _sigma_hadrons[0]->fill(_ecms);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        BinnedEstimatePtr<string> R;
        book(R,1,1,1+ix);
        divide(_sigma_hadrons[ix], _sigma_muons, R);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string>  _sigma_hadrons[2], _sigma_muons;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2006_I735496);

}
