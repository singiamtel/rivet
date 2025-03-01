#include "Rivet/Analysis.hh"
#include "Rivet/Tools/RivetPaths.hh"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace Rivet;
using namespace std;

class DummyAnalysis : public Analysis {
public:

  DummyAnalysis(string name = "DummyAnalysis")
    : Analysis(name)
  {	  }

  void analyze(const Event& event) { }

};

int main() {

  addAnalysisDataPath({"."});
  DummyAnalysis ana;
  ana.loadInfo();
  AnalysisInfo& anainfo = ana.info();
  assert(ana.name() == "DummyAnalysis");
  anainfo.setName("DummyAnalysisAlt");
  assert(ana.name() == "DummyAnalysisAlt");
  assert(ana.year() == "1234");
  assert(ana.experiment() == "DUMMYEX");
  assert(ana.collider() == "DUMMYCOLLIDER");
  assert(ana.inspireID() == "1234567");
  assert(ana.spiresID() == "98765");
  assert(ana.summary() == "Dummy info file");
  assert(ana.runInfo() == "Dummy run info");
  assert(ana.description() == "Dummy description.");
  assert(fuzzyEquals(ana.luminosity(), 123000.));
  assert(fuzzyEquals(ana.luminosityfb(), 123.0));
  vector<string> refs = ana.references();
  assert(refs.size() == 1);
  assert(refs[0] == "Dummy reference");
  assert(ana.bibKey() == "dummy");
  assert(ana.bibTeX() == "Dummy citation.");
  assert(ana.reentrant());
  assert(ana.status() == "VALIDATED");
  assert(ana.warning() == "dummy warning");
  vector<string> todos = ana.todos();
  assert(todos.size() == 1);
  assert(todos[0] == "dummy todo");
  vector<string> reltests = ana.validation();
  assert(reltests.size() == 1);
  assert(reltests[0] == "$A DUMMY");
  vector<string> keys = ana.keywords();
  assert(keys.size() == 1);
  assert(keys[0] == "key1");
  assert(ana.refMatch() == "x01");
  assert(ana.refUnmatch() == "y01");
  assert(ana.writerDoublePrecision() == "z01");
  assert(ana.refDataName() == "DummyAnalysisAlt");
  ana.setRefDataName("DummyRefName");
  assert(ana.refDataName() == "DummyRefName");
  ana.setRequiredBeamIDs({{2212,2212}});
  ana.setRequiredBeamEnergies({{7000.,7000.}});
  assert(ana.refFile() == "");
  vector<string> authors = ana.authors();
  assert(authors.size() == 1);
  assert(authors[0] == "Grogu <your@email.address>");
  ana.init();

  AnalysisInfo info2;
  assert(info2.name() == "");
  info2.setExperiment("DUMMY");
  info2.setYear("1234");
  assert(info2.name() == "");
  info2.setSpiresID("987654");
  assert(info2.name() == "DUMMY_1234_S987654");
  info2.setInspireID("123456");
  assert(info2.name() == "DUMMY_1234_I123456");

  return 0;
}
