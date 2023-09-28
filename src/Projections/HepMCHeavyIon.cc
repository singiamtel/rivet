// -*- C++ -*-
#include "Rivet/Projections/HepMCHeavyIon.hh"

#define IMPLEMENTATION(rettype, functionname, defret) \
rettype HepMCHeavyIon::functionname() const { \
  return _hi? _hi->functionname: defret;        \
}

namespace Rivet {


HepMCHeavyIon::HepMCHeavyIon() {
  setName("HepMCHeavyIon");
}

void HepMCHeavyIon::project(const Event& e) {
  _hi = e.genEvent()->heavy_ion();
  if ( !_hi )
    MSG_WARNING("Could not find the HepMC HeavyIon object");
}

IMPLEMENTATION(int, Ncoll_hard, -1)

IMPLEMENTATION(int, Npart_proj, -1)

IMPLEMENTATION(int, Npart_targ, -1)

IMPLEMENTATION(int, Ncoll, -1)

IMPLEMENTATION(int, N_Nwounded_collisions, -1)

IMPLEMENTATION(int, Nwounded_N_collisions, -1)

IMPLEMENTATION(int, Nwounded_Nwounded_collisions, -1)

IMPLEMENTATION(double, impact_parameter, -1.0)

IMPLEMENTATION(double, event_plane_angle, -1.0)

IMPLEMENTATION(double, sigma_inel_NN, -1.0)

IMPLEMENTATION(double, centrality, -1.0)

IMPLEMENTATION(double, user_cent_estimate, -1.0)

IMPLEMENTATION(int, Nspec_proj_n, -1)

IMPLEMENTATION(int, Nspec_targ_n, -1)

IMPLEMENTATION(int, Nspec_proj_p, -1)

IMPLEMENTATION(int, Nspec_targ_p, -1)

map<int,double> HepMCHeavyIon::participant_plane_angles() const {
  return _hi? _hi->participant_plane_angles: map<int,double>();
}

map<int,double> HepMCHeavyIon::eccentricities() const {
  return _hi? _hi->eccentricities: map<int,double>();
}

}
