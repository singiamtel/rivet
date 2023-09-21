// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Hadronic charged multiplicity measurement between 12 and 31.3 GeV
  class PLUTO_1980_I154270 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1980_I154270);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      for (const string label : {"9.4", "12.0", "13.0", "17.0", "22.0", "27.6", "30.2", "30.7", "31.3"}) {
        const double sqrts = std::stod(label);
        if (isCompatibleWithSqrtS(sqrts*GeV))  Ecm = label;
      }
      if (Ecm != "") {
        book(_mult, 1, 1, 1);
      }
      else {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV << " doesn't match any available analysis energy .");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      unsigned int nPart(0);
      for (const Particle& p : cfs.particles()) {
        // check if prompt or not
        ConstGenParticlePtr pmother = p.genParticle();
        ConstGenVertexPtr ivertex = pmother->production_vertex();
        bool prompt = true;
        while (ivertex) {
          vector<ConstGenParticlePtr> inparts = HepMCUtils::particles(ivertex, Relatives::PARENTS);
          int n_inparts = inparts.size();
          if (n_inparts < 1)  break;
          pmother = inparts[0]; // first mother particle
          int mother_pid = abs(pmother->pdg_id());
          if (mother_pid==PID::K0S || mother_pid==PID::LAMBDA) {
            prompt = false;
            break;
          }
          else if (mother_pid<6) {
            break;
          }
          ivertex = pmother->production_vertex();
        }
        if (prompt) ++nPart;
      }
      _mult->fill(Ecm, nPart);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_mult, 1./sumOfWeights());
    }

    /// @}


  private:

    BinnedProfilePtr<string> _mult;
    string Ecm = "";

  };


  RIVET_DECLARE_PLUGIN(PLUTO_1980_I154270);

}
