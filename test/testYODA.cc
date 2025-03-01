#include "Rivet/Tools/RivetYODA.hh"
#include "YODA/Histo.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Rivet;
using namespace std;


int main() {

  YODA::Histo1DPtr h1 = make_shared<YODA::Histo1D>(vector<double>{0.0, 0.5, 1.0});
  YODA::Histo1DPtr h2 = make_shared<YODA::Histo1D>(vector<double>{0.0, 0.5, 1.0});
  YODA::Histo1DPtr h3 = make_shared<YODA::Histo1D>(vector<double>{0.0, 1.0});
  YODA::Estimate1DPtr e1 = make_shared<YODA::Estimate1D>(vector<double>{0.0, 0.5, 1.0});

  assert( bookingCompatible(h1, h2) );
  assert( !bookingCompatible(h1, h3) );

  bool pass = false;
  try {
    copyAO<YODA::Histo1D>(h1, e1);
  }
  catch (const YODA::LogicError&) {
    pass = true;
  }
  if (!pass)  return EXIT_FAILURE;

  h1->setPath("histo");
  Multiplexer<YODA::Histo1D> hist({""}, *h1);
  assert( hist.basePath() == "/histo" );
  assert( hist.baseName() == "histo" );

  hist.setActiveWeightIdx(0);
  hist.newSubEvent();
  hist->fill(0.2, 0.25); // coords + weight fraction
  hist->fill(0.4, 0.25);
  hist.newSubEvent();
  hist->fill(0.75, 0.5);
  const vector<valarray<double>> weights = {{1.0}, {10.0}};
  hist.collapseEventGroup(weights, 0.1);
  assert( fuzzyEquals(hist.persistent(0)->bin(1).sumW(), 0.5) ); // 0.25 * 1.0 + 0.25 * 1.0
  assert( fuzzyEquals(hist.persistent(0)->bin(2).sumW(), 5.0) ); // 0.5 * 10.0

  Multiplexer<YODA::Counter> count({""}, YODA::Counter());

  count.setActiveWeightIdx(0);
  count.newSubEvent();
  count->fill(0.25); // weight fraction
  count.newSubEvent();
  count->fill(0.75);
  count.collapseEventGroup(weights, 0.1);
  assert( fuzzyEquals(count.persistent(0)->sumW(), 7.75) ); // 0.25 * 1 + 0.75 * 10

  return EXIT_SUCCESS;
}
