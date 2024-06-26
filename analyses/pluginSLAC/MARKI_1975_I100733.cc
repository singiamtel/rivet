// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement 
  class MARKI_1975_I100733 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1975_I100733);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(FinalState(), "FS");
      // Book histograms
      for(unsigned int ix=0;ix<2;++ix) {
        if(ix==0) {
          book(_c_hadrons[ix], 1+ix,1,1);
        }
        else {
          book(_c_hadrons[ix], "/TMP/sigma_hadrons_"+toString(ix), refData<YODA::BinnedEstimate<string>>(1+ix,1,1));
          book(_c_muons,   "/TMP/sigma_muons_"  +toString(ix), refData<YODA::BinnedEstimate<string>>(1+ix,1,1));
        }
        for (const string& en : _c_hadrons[ix].binning().edges<0>()) {
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
      if (isCompatibleWithSqrtS(3*GeV)) {
      	book(_h_charged, 3, 1, 1);
      }
      else if (isCompatibleWithSqrtS(4.8*GeV)) {
        book(_h_charged, 3, 1, 2);
      }
      else if (isCompatibleWithSqrtS(7.4*GeV)) {
        book(_h_charged, 3, 1, 3);
      }
      if(_ecms[0].empty() &&_ecms[1].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
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
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
      	_c_muons->fill(_ecms[1]); // mu+mu- + photons
      }
      // everything else
      else {
	for(unsigned int ix=0;ix<2;++ix) _c_hadrons[ix]->fill(_ecms[ix]);
        if(_h_charged) {
          if(_edges.empty()) _edges = _h_charged->xEdges();
          for (const Particle& p : fs.particles()) {
            if(PID::isCharged(p.pid())) {
              double x = 2.*p.p3().mod()/sqrtS();
              const size_t idx = _axis.index(x);
              if(idx && idx <=_edges.size())
                _h_charged->fill(_edges[idx-1]);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_h_charged) {
      	scale(_h_charged, crossSection()/ sumOfWeights() /microbarn*sqr(sqrtS()));
        for(auto & b : _h_charged->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axis.width(idx));
        }
      }
      const double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale(_c_hadrons[0], fact);
      BinnedEstimatePtr<string> mult;
      book(mult,2,1,1);
      divide(_c_hadrons[1],_c_muons,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons[2], _c_muons;
    string _ecms[2];
    BinnedHistoPtr<string> _h_charged;
    vector<string> _edges;
    YODA::Axis<double> _axis{0.,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,
        0.32,0.34,0.36,0.38,0.40,0.44,0.48,0.52,0.56,0.64,0.72,0.80,0.88};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKI_1975_I100733);


}
