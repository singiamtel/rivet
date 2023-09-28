// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @author Andy Buckley
  class MC_PRINTEVENT : public Analysis {
  public:

    /// Constructor
    MC_PRINTEVENT()
      : Analysis("MC_PRINTEVENT")
    {    }


    /// @name Analysis methods
    /// @{

    void init() {

      /// Set up particle name map
      // quarks & gluons
      _pnames[1] = "d";
      _pnames[-1] = "d~";
      _pnames[2] = "u";
      _pnames[-2] = "u~";
      _pnames[3] = "s";
      _pnames[-3] = "s~";
      _pnames[4] = "c";
      _pnames[-4] = "c~";
      _pnames[5] = "b";
      _pnames[-5] = "b~";
      _pnames[6] = "t";
      _pnames[-6] = "t~";
      // bosons
      _pnames[21] = "g";
      _pnames[22] = "gamma";
      _pnames[23] = "Z0";
      _pnames[24] = "W+";
      _pnames[-24] = "W-";
      _pnames[25] = "h0";
      _pnames[26] = "h0";
      // leptons
      _pnames[11] = "e-";
      _pnames[-11] = "e+";
      _pnames[13] = "mu-";
      _pnames[-13] = "mu+";
      _pnames[15] = "tau-";
      _pnames[-15] = "tau+";
      _pnames[12] = "nu_e";
      _pnames[-12] = "nu_e~";
      _pnames[14] = "nu_mu";
      _pnames[-14] = "nu_mu~";
      _pnames[16] = "nu_tau";
      _pnames[-16] = "nu_tau~";
      // common hadrons
      _pnames[111] = "pi0";
      _pnames[211] = "pi+";
      _pnames[-211] = "pi-";
      _pnames[221] = "eta";
      _pnames[331] = "eta'";
      _pnames[113] = "rho0";
      _pnames[213] = "rho+";
      _pnames[-213] = "rho-";
      _pnames[223] = "omega";
      _pnames[333] = "phi";
      _pnames[130] = "K0L";
      _pnames[310] = "K0S";
      _pnames[311] = "K0";
      _pnames[-311] = "K0";
      _pnames[321] = "K+";
      _pnames[-321] = "K-";
      _pnames[313] = "K*0";
      _pnames[-313] = "K*0~";
      _pnames[323] = "K*+";
      _pnames[-323] = "K*-";
      _pnames[411] = "D+";
      _pnames[-411] = "D-";
      _pnames[421] = "D0";
      _pnames[-421] = "D0~";
      _pnames[413] = "D*+";
      _pnames[-413] = "D*-";
      _pnames[423] = "D*0";
      _pnames[-423] = "D*0~";
      _pnames[431] = "Ds+";
      _pnames[-431] = "Ds-";
      _pnames[433] = "Ds*+";
      _pnames[-433] = "Ds*-";
      _pnames[511] = "B0";
      _pnames[-511] = "B0~";
      _pnames[521] = "B+";
      _pnames[-521] = "B-";
      _pnames[513] = "B*0";
      _pnames[-513] = "B*0~";
      _pnames[523] = "B*+";
      _pnames[-523] = "B*-";
      _pnames[531] = "B0s";
      _pnames[541] = "Bc+";
      _pnames[-541] = "Bc-";
      _pnames[441] = "eta_c(1S)";
      _pnames[443] = "J/psi(1S)";
      _pnames[551] = "eta_b(1S)";
      _pnames[553] = "Upsilon(1S)";
      _pnames[2212] = "p+";
      _pnames[-2212] = "p-";
      _pnames[2112] = "n";
      _pnames[-2112] = "n~";
      _pnames[2224] = "Delta++";
      _pnames[2214] = "Delta+";
      _pnames[2114] = "Delta0";
      _pnames[1114] = "Delta-";
      _pnames[3122] = "Lambda";
      _pnames[-3122] = "Lambda~";
      _pnames[3222] = "Sigma+";
      _pnames[-3222] = "Sigma+~";
      _pnames[3212] = "Sigma0";
      _pnames[-3212] = "Sigma0~";
      _pnames[3112] = "Sigma-";
      _pnames[-3112] = "Sigma-~";
      _pnames[4122] = "Lambda_c+";
      _pnames[-4122] = "Lambda_c-";
      _pnames[5122] = "Lambda_b";
      // exotic
      _pnames[32] = "Z'";
      _pnames[34] = "W'+";
      _pnames[-34] = "W'-";
      _pnames[35] = "H0";
      _pnames[36] = "A0";
      _pnames[37] = "H+";
      _pnames[-37] = "H-";
      // shower-specific
      _pnames[91] = "cluster";
      _pnames[92] = "string";
      _pnames[9922212] = "remn";
      _pnames[1103] = "dd";
      _pnames[2101] = "ud0";
      _pnames[2103] = "ud1";
      _pnames[2203] = "uu";

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      RivetHepMC::Print::content(*(event.genEvent()));

    }


    /// Normalise histograms etc., after the run
    void finalize() {    }

    /// @}

  private:

    map<long, string> _pnames;


  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_PRINTEVENT);

}
