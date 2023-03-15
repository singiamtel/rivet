#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/TypeRegistry.hh"
#include "YODA/BinnedDbn.h"
#include "YODA/IO.h"

// #include <regex>
#include <sstream>
using namespace std;

namespace Rivet {

  string getDatafilePath(const string& papername) {
    /// Try to find a YODA file matching this analysis name
    const string path1 = findAnalysisRefFile(papername + ".yoda");
    if (!path1.empty()) return path1;
    const string path2 = findAnalysisRefFile(papername + ".yoda.gz");
    if (!path2.empty()) return path2;
    throw Rivet::Error("Couldn't find a ref data file for '" + papername +
                       "' in data path, '" + getRivetDataPath() + "', or '.'");
  }


  map<string, YODA::AnalysisObjectPtr> getRefData(const string& papername) {
    const string datafile = getDatafilePath(papername);

    // Read the data objects
    vector<YODA::AnalysisObject*> aovec;
    YODA::read(datafile, aovec);

    // Return value, to be populated
    map<string, YODA::AnalysisObjectPtr> rtn;
    for (YODA::AnalysisObject* ao : aovec) {
      YODA::AnalysisObjectPtr refdata(ao);
      if (!refdata) continue;
      const string plotpath = refdata->path();
      // Split path at "/" and only return the last field, i.e. the histogram ID
      const size_t slashpos = plotpath.rfind("/");
      const string plotname = (slashpos+1 < plotpath.size()) ? plotpath.substr(slashpos+1) : "";
      rtn[plotname] = refdata;
    }
    return rtn;
  }


  /// @todo Switch to this path-handling solution based on std::regex
  //
  // AOPath::AOPath(string fullpath) {
  // // First check if this is a global system path
  // _path = fullpath;
  // std::regex resys("^(/RAW)?/([^\\[/]+)(\\[(.+)\\])?$");
  // smatch m;
  // _valid = regex_search(fullpath, m, resys);
  // if (_valid) {
  // _raw = (m[1] == "/RAW");
  // _name = m[2];
  // _weight = m[4];
  // return;
  // }
  // // If not, assume it is a normal analysis path.
  // std::regex repath("^(/RAW)?(/REF)?/([^/:]+)(:[^/]+)?(/TMP)?/([^\\[]+)(\\[(.+)\\])?");
  // _valid = regex_search(fullpath, m, repath);
  // if ( !_valid ) return;
  // _raw = (m[1] == "/RAW");
  // _ref = (m[2] == "/REF");
  // _analysis = m[3];
  // _optionstring = m[4];
  // _tmp = (m[5] == "/TMP");
  // _name = m[6];
  // _weight = m[8];
  // std::regex reopt(":([^=]+)=([^:]+)");
  // string s = _optionstring;
  // while ( regex_search(s, m, reopt) ) {
  // _options[m[1]] = m[2];
  // s = m.suffix();
  // }
  // }


  bool AOPath::init(string fullpath) {
    if ( fullpath.substr(0,5) == "/RAW/" ) {
      _raw = true;
      return init(fullpath.substr(4));
    }
    if ( fullpath.substr(0,5) == "/REF/" ) {
      _ref = true;
      return init(fullpath.substr(4));
    }
    if ( fullpath[0] != '/' ) return false;

    fullpath = fullpath.substr(1);

    if ( fullpath.size() < 2 ) return false;

    if ( !chopweight(fullpath) ) return false;

    string::size_type p = fullpath.find("/");
    if ( p == 0 ) return false;
    if ( p == string::npos ) {
      _name = fullpath;
      return true;
    }
    _analysis = fullpath.substr(0, p);
    _name = fullpath.substr(p + 1);

    if ( _name.substr(0, 4) == "TMP/" ) {
      _name = _name.substr(4);
      _tmp = true;
    }

    if ( !chopoptions(_analysis) ) return false;

    fixOptionString();

    return true;
  }

  bool AOPath::chopweight(string & fullpath) {
    if ( fullpath.back() != ']' ) return true;
    string::size_type p = fullpath.rfind("[");
    if ( p == string::npos ) return false;
    _weight = fullpath.substr(p + 1);
    _weight.pop_back();
    fullpath = fullpath.substr(0, p);
    return true;
  }

  bool AOPath::chopoptions(string & anal) {
    string::size_type p = anal.rfind(":");
    if ( p == string::npos ) return true;
    string opts = anal.substr(p + 1);
    string::size_type pp = opts.find("=");
    if ( pp == string::npos ) return false;
    _options[opts.substr(0, pp)] = opts.substr(pp + 1);
    anal = anal.substr(0, p);
    return chopoptions(anal);
  }

  void AOPath::fixOptionString() {
    ostringstream oss;
    for ( auto optval : _options )
      oss << ":" << optval.first << "=" << optval.second;
    _optionstring = oss.str();
  }

  string AOPath::mkPath() const {
    ostringstream oss;
    if ( isRaw() ) oss << "/RAW";
    else if ( isRef() ) oss << "/REF";
    if ( _analysis != "" ) oss << "/" << analysis();
    for ( auto optval : _options )
      oss << ":" << optval.first << "=" << optval.second;
    if ( isTmp() ) oss << "/TMP";
    oss << "/" << name();
    if ( weight() != "" )
      oss << "[" << weight() << "]";
    return oss.str();
  }

  void AOPath::debug() const {
    cout << "Full path:  " << _path << endl;
    if ( !_valid ) {
      cout << "This is not a valid analysis object path" << endl << endl;
      return;
    }
    cout << "Check path: " << mkPath() << endl;
    cout << "Analysis:   " << _analysis << endl;
    cout << "Name:       " << _name << endl;
    cout << "Weight:     " << _weight << endl;
    cout << "Properties: ";
    if ( _raw ) cout << "raw ";
    if ( _tmp ) cout << "tmp ";
    if ( _ref ) cout << "ref ";
    cout << endl;
    cout << "Options:    ";
    for ( auto opt : _options )
      cout << opt.first << "->" << opt.second << " ";
    cout << endl << endl;
  }

}
