// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALEPH_2014_I1267648 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_2014_I1267648);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=0;ix<5;++ix)
        book(_h[ix], "TMP/h_"+toString(ix+1),refData(1+ix,1,1));
      book(_c,"TMP/ntau");
    }


    void findDecayProducts(const Particle &mother, const int isign, unsigned int &nstable, unsigned int &npip,
                           unsigned int &npim, unsigned int &npi0, FourMomentum &ptot) {
      for (const Particle &p : mother.children()) {
        int id = p.pid();
        if (id == PID::KPLUS || id == PID::KMINUS) {
          ++nstable;
          ptot += p.momentum();
        }
        else if (id*isign == PID::PIPLUS) {
          ++npip;
          ++nstable;
          ptot += p.momentum();
        }
        else if (id*isign == PID::PIMINUS) {
          ++npim;
          ++nstable;
          ptot += p.momentum();
        }
        else if (id == PID::PI0) {
          ++nstable;
          ++npi0;
          ptot += p.momentum();
        }
        else if (id == PID::PHOTON)  continue;
        else if (!p.children().empty())  findDecayProducts(p, isign, nstable, npip, npim, npi0, ptot);
        else  ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for (const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        _c->fill();
        FourMomentum ptot;
        unsigned int nstable(0), npip(0), npim(0), npi0(0);
        findDecayProducts(tau,tau.pid()/tau.abspid(),nstable,npip,npim,npi0,ptot);
        // tau -> pi pi0 nu_tau (both charges)
        if (npim==1 && npi0==1 && nstable==3)  _h[0]->fill(ptot.mass2());
        // tau -> pi pi0 pi0 nu_tau (both charges)
        else if (npim==1 && npi0==2 && nstable==4)  _h[1]->fill(ptot.mass2());
        //    tau -> pi pi0 pi0 pi0         (3,1,1)
        else if (npim==1 && npi0==3 && nstable==5)  _h[2]->fill(ptot.mass2());
        //    tau -> 3 charged pions        (4,1,1)
        else if (npim==2 && npip==1 && nstable==4)  _h[3]->fill(ptot.mass2());
        //    tau -> 3 charged pions + pi0  (5,1,1)
        else if (npim==2 && npip==1 && npi0==1 && nstable==5)  _h[4]->fill(ptot.mass2());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<5;++ix) {
        scale(_h[ix] ,100./ *_c);
        Estimate1DPtr tmp;
        book(tmp,1+ix,1,1);
        barchart(_h[ix],tmp);
      }
    }

    /// @}


  private:


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    CounterPtr _c;
    /// @}
  };


  RIVET_DECLARE_PLUGIN(ALEPH_2014_I1267648);

}
