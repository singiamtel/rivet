// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Cross sections for $e^+e^-\to J/\psi + \pi^+\pi^-$, $\pi^0\pi^0$ or $K^+K^-$
  class CLEO_2006_I710864 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2006_I710864);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      _mode = 0;
      if ( getOption("MODE") == "SIGMA" ) _mode = 0;
      else if ( getOption("MODE") == "DECAY" ) _mode = 1;
      // set the PDG code
      _pid = getOption<int>("PID", 9030443);
      // setup for cross sections
      if (_mode==0) {
        // projections
        declare(FinalState(), "FS");
        declare(UnstableParticles(Cuts::pid== 443), "UFS");
        // histograms
        for (unsigned int ix=0; ix<3; ++ix) {
          book(_sigma[ix],1,1,1+ix);
        }
      for (const string& en : _sigma[0].binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms = en;
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms = en;
            break;
          }
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
      // setup for decays
      else {
        // projections
        UnstableParticles ufs = UnstableParticles(Cuts::abspid==_pid);
        declare(ufs, "UFS");
        DecayedParticles PSI(ufs);
        PSI.addStable(443);
        declare(PSI, "PSI");
        book(_h,2,1,1);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle& child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // cross section
      if (_mode==0) {
        const FinalState& fs = apply<FinalState>(event, "FS");
        map<long,int> nCount;
        int ntotal(0);
        for (const Particle& p : fs.particles()) {
          nCount[p.pid()] += 1;
          ++ntotal;
        }
        for (const Particle& psi : apply<UnstableParticles>(event, "UFS").particles()) {
          if (psi.children().empty()) continue;
          map<long,int> nRes = nCount;
          int ncount = ntotal;
          findChildren(psi,nRes,ncount);
          if (ncount!=2) continue;
          if (nRes.find(211)!=nRes.end() && nRes.find(-211)!=nRes.end() &&
             nRes[211]==1 && nRes[-211]==1) {
            _sigma[0]->fill(_ecms);
            break;
          }
          else if (nRes.find(111)!=nRes.end() && nRes[111]==2) {
            _sigma[1]->fill(_ecms);
            break;
          }
          else if (nRes.find(321)!=nRes.end() && nRes.find(-321)!=nRes.end() &&
                   nRes[321]==1 && nRes[-321]==1) {
            _sigma[2]->fill(_ecms);
            break;
          }
        }
      }
      // decays
      else {
        DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
        for (unsigned int ix=0; ix<PSI.decaying().size(); ++ix) {
          if (!PSI.modeMatches(ix,3,mode)) continue;
          const Particle & pip  = PSI.decayProducts()[ix].at( 211)[0];
          const Particle & pim  = PSI.decayProducts()[ix].at(-211)[0];
          double mpipi = (pip.mom()+pim.mom()).mass();
          _h->fill(mpipi);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // sigma
      if (_mode==0) {
        scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
      }
      // decay
      else {
        normalize(_h, 1.0, false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    int _mode;
    int _pid;
    Histo1DPtr _h;
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    const map<PdgId,unsigned int> mode = { { 211,1}, {-211,1}, { 443,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2006_I710864);

}
