#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include <fstream>

using namespace std;

int main() {

  Rivet::AnalysisHandler _ah[2];
  Rivet::Log::setLevel("Rivet", Rivet::Log::DEBUG);
  Rivet::addAnalysisLibPath(".");

  _ah[0].addAnalyses({{ "MC_JETS"s }});
  _ah[1].addAnalyses({{ "MC_JETS"s }});

  std::vector<std::string> outnames{ "out-merge0.yoda"s, "out-merge1.yoda"s };

  auto reader = Rivet::RivetHepMC::deduce_reader("testAPI.hepmc");
  std::shared_ptr<Rivet::GenEvent> evt = make_shared<Rivet::GenEvent>();
  size_t numEvents = 0;
  while ( Rivet::HepMCUtils::readEvent(reader, evt) ) {
    ++numEvents;
    // Analyse current event
    _ah[ int(numEvents > 5) ].analyze(*evt);
  }

  double xs, xserr;
  vector<vector<double>> data;
  for (size_t i = 0; i < 2; ++i) {
    data.push_back(_ah[i].serializeContent(true));
    _ah[i].finalize();
    _ah[i].writeData(outnames[i]);
    if (i) {
      xs = _ah[i].nominalCrossSection();
      xserr = _ah[i].nominalCrossSectionError();
    }
  }

  // perform a dummy MPI all-reduce
  vector<double> d; d.reserve(data[0].size());
  for (size_t i = 0; i < data[0].size(); ++i) {
    d.push_back(data[0][i] + data[1][i]);
  }

  // Test AH merging
  _ah[0].merge(_ah[1]);
  _ah[0].finalize();
  _ah[0].writeData("out-merge-ahmerge.yoda");

  // Test deserialisation
  _ah[0].deserializeContent(d, 2);
  _ah[0].setCrossSection(xs, xserr);
  _ah[0].finalize();
  _ah[0].writeData("out-merge-serial.yoda");

  // Test equivalent merging
  Rivet::AnalysisHandler _ahmerge[2];
  _ahmerge[0].mergeYODAs(outnames, {}, {}, {}, {}, true);
  _ahmerge[0].finalize();
  _ahmerge[0].writeData("out-merge-equiv.yoda");

  // Test non-equivalent merging
  _ahmerge[1].mergeYODAs(outnames, {}, {}, {}, {}, false);
  _ahmerge[1].finalize();
  _ahmerge[1].writeData("out-merge-nonequiv.yoda");

  return 0;
}

