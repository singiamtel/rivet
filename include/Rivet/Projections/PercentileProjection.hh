// -*- C++ -*-
#ifndef RIVET_PERCENTILEPROJECTION_HH
#define RIVET_PERCENTILEPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <map>

namespace Rivet {

  enum class PercentileOrder { INCREASING, DECREASING };
  

  /// @brief class for projections that reports the percentile for a
  /// given SingleValueProjection when initialized with a Histo1D of the
  /// distribution in the SingleValueProjection.
  ///
  /// @author Leif Lönnblad
  class PercentileProjection : public SingleValueProjection {
  public:

    using SingleValueProjection::operator =;

    /// Constructor taking a SingleValueProjection and a calibration
    /// histogram. If increasing it means that low values corresponds to
    /// lower percentiles.
    ///
    /// @todo Use mkScatter to pass this to the Scatter2D-calibrated version?
    PercentileProjection(const SingleValueProjection & sv, const Histo1D& calhist,
                         PercentileOrder pctorder=PercentileOrder::DECREASING)
      : _calhist("EMPTY"),
	_increasing(pctorder == PercentileOrder::INCREASING)
    {
      setName("PercentileProjection");
      declare(sv, "OBSERVABLE");
      //if ( !calhist ) return;
      MSG_DEBUG("Constructing PercentileProjection from " << calhist.path());
      _calhist = calhist.path();
      int N = calhist.numBins();
      double sum = calhist.sumW();
      if ( _increasing ) {
        double acc = 0.0;
        for (int i = 0; i <= N; ++i) {
          acc += calhist.bin(i).sumW();
          _table.insert(make_pair(calhist.bin(i).xMax(), 100.0*acc/sum));
        }
      }
      else {
        double acc = 0.0;
        for (int i = N+1; i > 0; --i) {
          acc += calhist.bin(i).sumW();
          _table.insert(make_pair(calhist.bin(i).xMin(), 100.0*acc/sum));
        }
      }
      if (getLog().isActive(Log::DEBUG)) {
        MSG_DEBUG("Mapping from observable to percentile:");
        for (auto p : _table) {
          std::cout << std::setw(16) << p.first << " -> "
                    << std::setw(16) << p.second << "%" << std::endl;
          if (not _increasing and p.second <= 0) break;
          if (_increasing and p.second >= 100) break;
        }
      }
    }

    
    // Constructor taking a SingleValueProjection and a calibration
    // histogram. If increasing it means that low values corresponds to
    // lower percentiles.
    PercentileProjection(const SingleValueProjection & sv, const Estimate1D& calest,
                         PercentileOrder pctorder=PercentileOrder::DECREASING)
      : _calhist("EMPTY"),
	_increasing(pctorder == PercentileOrder::INCREASING)
    {
      declare(sv, "OBSERVABLE");

      //if ( !calest ) return;
      MSG_DEBUG("Constructing PercentileProjection from " << calest.path());
      _calhist = calest.path();
      int N = calest.numBins();
      double sum = 0.0;
      for (const auto &b : calest.bins() )  sum += b.val();

      double acc = 0.0;
      if ( _increasing ) {
        _table.insert(make_pair(calest.bin(1).xMin(), 100.0*acc/sum));
        for ( int i = 0; i < N; ++i ) {
          acc += calest.bin(i+1).val();
          _table.insert(make_pair(calest.bin(i+1).xMax(), 100.0*acc/sum));
        }
      }
      else {
        _table.insert(make_pair(calest.bin(N).xMax(), 100.0*acc/sum));
        for (int i = N - 1; i >= 0; --i) {
          acc += calest.bin(i+1).val();
          _table.insert(make_pair(calest.bin(i+1).xMin(), 100.0*acc/sum));
        }
      }
    }

    
    RIVET_DEFAULT_PROJ_CLONE(PercentileProjection);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    // The projection function takes the assigned SingeValueProjection
    // and sets the value of this projection to the corresponding
    // percentile. If no calibration has been provided, -1 will be
    // returned. If values are outside of the calibration histogram, 0
    // or 100 will be returned.
    void project(const Event& e) {
      clear();
      if ( _table.empty() ) return;
      auto&  pobs = apply<SingleValueProjection>(e, "OBSERVABLE");
      double obs  = pobs();
      double pcnt = lookup(obs);
      if ( pcnt >= 0.0 ) set(pcnt);
      MSG_DEBUG("Observable(" << pobs.name() << ")="
                << std::setw(16) << obs
                << "-> Percentile=" << std::setw(16) << pcnt << "%");
    }

    // Standard comparison function.
    CmpState compare(const Projection& p) const {
      const PercentileProjection pp = dynamic_cast<const PercentileProjection&>(p);
      return mkNamedPCmp(p, "OBSERVABLE") ||
        cmp(_increasing, pp._increasing) ||
        cmp(_calhist, pp._calhist);
    }


  protected:

    // The (interpolated) lookup table
    double lookup(double obs) const {
      auto low = _table.upper_bound(obs);
      if ( low == _table.end() ) return _increasing? 100.0: 0.0;
      if ( low == _table.begin() ) return _increasing? 0.0: 100.0;
      auto high = low--;
      return low->second + (obs - low->first)*(high->second - low->second)/
        (high->first - low->first);
    }

    // Astring identifying the calibration histogram.
    string _calhist;

    // A lookup table to find (by interpolation) the percentile given
    // the value of the underlying SingleValueProjection.
    map<double,double> _table;

    // A flag to say whether the distribution should be integrated from
    // below or above.
    bool _increasing;

  };


}

#endif
