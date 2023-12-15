// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2016_I1389855 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1389855);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_nhb1, 1, 1, 3);
      book(_nhb2, 1, 1, 4);

      // Check which indices match the Ecm energy
      for (const auto& est : refData<YODA::BinnedEstimate<int>>(1, 1, 1).bins()) {
        if (inRange(sqrtS()/MeV, est.valMin(), est.valMax())) {
          edges.push_back(est.xEdge());
        }
      }
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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

      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if (p.children().empty()) continue;
        // find the omega
        if (p.pid()== 10553|| p.pid()==110553) {
          map<long,int> nRes = nCount;
          int ncount = ntotal;
          findChildren(p,nRes,ncount);
          // omega pi+pi-
          if(ncount!=2) continue;
          bool matched = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==211) {
              if (val.second !=1) {
                matched = false;
                break;
              }
            }
            else if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            for (const int Ecm : edges) {
              if (p.pid() == 10553) {
                _nhb1->fill(Ecm);
              }
              else {
                _nhb2->fill(Ecm);
              }
            }
            break;
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/sumOfWeights()/picobarn;
      scale({_nhb1, _nhb2}, sf);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _nhb1, _nhb2;
    vector<int> edges;
    /// @}

  };

  RIVET_DECLARE_PLUGIN(BELLE_2016_I1389855);
}
