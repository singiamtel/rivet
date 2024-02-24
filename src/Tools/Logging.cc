#include "Rivet/Tools/Logging.hh"
#include <ctime>
#include <unistd.h>
using namespace std;

namespace Rivet {


  thread_local Log::LogMap Log::existingLogs;
  thread_local Log::LevelMap Log::defaultLevels;
  bool Log::showTimestamp = false;
  bool Log::showLogLevel = true;
  bool Log::showLoggerName = true;
  bool Log::useShellColors = true;
  const int Log::END_COLOR;


  Log::Log(const string& name)
    : _name(name), _level(INFO) { }


  Log::Log(const string& name, int level)
    : _name(name), _level(level) { }


  /// @todo Add single static setLevel
  void _updateLevels(const Log::LevelMap& defaultLevels, Log::LogMap& existingLogs) {
    /// @todo Check ordering - "Foo" should come before "Foo.Bar"
    for (Log::LevelMap::const_iterator lev = defaultLevels.begin(); lev != defaultLevels.end(); ++lev) {
      for (Log::LogMap::iterator log = existingLogs.begin(); log != existingLogs.end(); ++log) {
        if (log->first.find(lev->first) == 0) {
          log->second.setLevel(lev->second);
        }
      }
    }
  }


  void Log::setLevel(const string& name, int level) {
    defaultLevels[name] = level;
    //cout << name << " -> " << level << '\n';
    _updateLevels(defaultLevels, existingLogs);
  }


  void Log::setLevels(const LevelMap& logLevels) {
    for (LevelMap::const_iterator lev = logLevels.begin(); lev != logLevels.end(); ++lev) {
      defaultLevels[lev->first] = lev->second;
    }
    _updateLevels(defaultLevels, existingLogs);
  }


  Log& Log::getLog(const string& name) {
    auto theLog = existingLogs.find(name);
    if (theLog == existingLogs.end()) {
      int level = INFO;
      // Try running through all parent classes to find an existing level
      string tmpname = name;
      bool triedAllParents = false;
      while (! triedAllParents) {
        // Is there a default level?
        if (defaultLevels.find(tmpname) != defaultLevels.end()) {
          level = defaultLevels.find(tmpname)->second;
          break;
        }
        // Is there already such a logger? (NB. tmpname != name in later iterations)
        if (existingLogs.find(tmpname) != existingLogs.end()) {
          level = existingLogs.find(tmpname)->second.getLevel();
          break;
        }
        // Crop the string back to the next parent level
        size_t lastDot = tmpname.find_last_of(".");
        if (lastDot != string::npos) {
          tmpname = tmpname.substr(0, lastDot);
        } else {
          triedAllParents = true;
        }
      }
      // for (LevelMap::const_iterator l = defaultLevels.begin(); l != defaultLevels.end(); ++l) {
      //
      // }

      // emplace returns pair<iterator,bool>
      auto result = existingLogs.emplace(name, Log(name, level));
      theLog = result.first;
    }
    return theLog->second;
  }


  string Log::getLevelName(int level) {
    switch(level) {
    case TRACE:
      return "TRACE";
    case DEBUG:
      return "DEBUG";
    case INFO:
      return "INFO";
    case WARN:
      return "WARN";
    case ERROR:
      return "ERROR";
    case CRITICAL:
      return "CRITICAL";
    default:
      return "";
    }
  }


  string Log::getColorCode(int level) {
    // Skip codes if
    if (!Log::useShellColors) return "";
    const static bool IS_TTY = isatty(1);
    if (!IS_TTY) return "";

    static const ColorCodes TTY_CODES = {
      {TRACE, "\033[0;36m"},
      {DEBUG, "\033[0;34m"},
      {INFO, "\033[0;32m"},
      {WARN, "\033[0;33m"},
      {ERROR, "\033[0;31m"},
      {CRITICAL, "\033[0;31m"},
      {END_COLOR, "\033[0m"} // end-color code
    };
    try {
      return TTY_CODES.at(level);
    } catch (...) {
      return "";
    }
  }


  Log::Level Log::getLevelFromName(const string& level) {
    if (level == "TRACE") return TRACE;
    if (level == "DEBUG") return DEBUG;
    if (level == "INFO") return INFO;
    if (level == "WARN") return WARN;
    if (level == "ERROR") return ERROR;
    if (level == "CRITICAL") return CRITICAL;
    throw Error("Couldn't create a log level from string '" + level + "'");
  }


  string Log::formatMessage(int level, const string& message) {
    string out;
    out += getColorCode(level);

    if (Log::showLoggerName) {
      out += getName();
      out += ": ";
    }

    if (Log::showLogLevel) {
      out += Log::getLevelName(level);
      out += " ";
    }

    if (Log::showTimestamp) {
      time_t rawtime;
      time(&rawtime);
      char* timestr = ctime(&rawtime);
      timestr[24] = ' ';
      out += timestr;
      out += " ";
    }

    out += getColorCode(END_COLOR);
    out += " ";
    out += message;

    return out;
  }


  void Log::log(int level, const string& message) {
    if (isActive(level)) {
        if (level > WARNING)
	      cerr << formatMessage(level, message) << '\n';
  	else
	      cout << formatMessage(level, message) << '\n';
    }
  }


  ostream& operator<<(Log& log, int level) {
    if (log.isActive(level)) {
    	if (level > Log::WARNING) {
      		cerr << log.formatMessage(level, "");
      		return cerr;
      	} else {
      		cout << log.formatMessage(level, "");
      		return cout;
      	}
    } else {
      static ostream devNull(nullptr);
      return devNull;
    }
  }

}
