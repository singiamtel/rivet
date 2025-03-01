#include "Rivet/Tools/HistoGroup.hh"
//#include "YODA/Histo.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Rivet;
using namespace std;


int main() {

  HistoGroup<double,double> group(vector<double>{0.0, 0.5, 1.0});

  for (auto& b : group.bins()) {
    Histo1DPtr hist({""}, YODA::Histo1D(vector<double>{0.0, 1.0}));
    b = hist;
    b.get()->newSubEvent();
  }

  group.fill(0.25, 0.5);

  for (auto& b : group.bins()) {
    b.get()->collapseEventGroup(vector<valarray<double>>{{1.0}});
    b.get()->newSubEvent();
  }

  group.fill(0.75, 0.5);

  for (auto& b : group.bins()) {
    b.get()->collapseEventGroup(vector<valarray<double>>{{10.0}});
    b.get()->setActiveWeightIdx(0);
  }

  assert(group.fillDim() == 2);

  assert( fuzzyEquals(group.numEntries(), 2.0) );
  assert( fuzzyEquals(group.effNumEntries(), 2.0) );
  assert( fuzzyEquals(group.sumW(), 11.0) );
  assert( fuzzyEquals(group.sumW2(), 101.0) );

  group.scaleW(2.0);
  assert( fuzzyEquals(group.sumW(), 22.0) );

  group.reset();

  assert( group.numBins() == 2 );
  for (auto& b : group.bins()) {
    assert( b.get() == nullptr );
  }

  return EXIT_SUCCESS;
}
