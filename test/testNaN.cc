#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


class NanTest : public Rivet::Analysis {
public:

  RIVET_DEFAULT_ANALYSIS_CTOR(NanTest);

  void init() {
    book(_h_test, "test", 50, 66.0, 116.0);
  }

  void analyze(const Rivet::Event& /* e */) {
    cout << "Normal fill" << endl;
    _h_test->fill(90.);

    cout << "Underflow fill" << endl;
    _h_test->fill(30.);

    cout << "Overflow fill" << endl;
    _h_test->fill(130.);

    cout << "-Inf fill" << endl;
    int res = _h_test->fill(-numeric_limits<double>::infinity());
    if (res != 0) throw;

    cout << "Inf fill" << endl;
    res = _h_test->fill(numeric_limits<double>::infinity());
    if (res != (int)_h_test->numBins(true)-1) throw;

    cout << "NaN fill" << endl;
    res = _h_test->fill(numeric_limits<double>::quiet_NaN());
    if (res != -1) throw;
  }

private:
  Rivet::Histo1DPtr _h_test;
};

RIVET_DECLARE_PLUGIN(NanTest);

int main() {

  Rivet::AnalysisHandler rivet;
  rivet.addAnalysis("NanTest");

  auto reader = Rivet::RivetHepMC::deduce_reader("testAPI.hepmc");
  std::shared_ptr<Rivet::GenEvent> evt = make_shared<Rivet::GenEvent>();
  double sum_of_weights = 0.0;

  while ( Rivet::HepMCUtils::readEvent(reader, evt) ) {
    // Analyse current event
    rivet.analyze(*evt);
    sum_of_weights += evt->weights()[0];
  }
  (void)sum_of_weights; // suppress unused variable warning

  rivet.setCrossSection(make_pair(1.0, 0.1));
  rivet.finalize();
  rivet.writeData("NaN.yoda");

  return 0;
}
