#ifndef RIVET_Cutflow_HH
#define RIVET_Cutflow_HH

#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Utils.hh"

#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

namespace Rivet {

  namespace {
    template<typename T>
    std::vector<T> operator+(std::vector<T>&& res, const std::vector<T>& vadd) {
      res.insert(std::end(res), std::begin(vadd), std::end(vadd));
      return std::move(res);
    }
  }

  /// A tracker of numbers & fractions of events passing sequential cuts
  class Cutflow : public YODA::BinnedHisto<std::string> {
  public:

    using BaseT = YODA::BinnedHisto<std::string>;
    using BinningT = BaseT::BinningT;
    using FillType = BaseT::FillType;
    using BinType = BaseT::BinT;
    using Ptr = std::shared_ptr<Cutflow>;
    using AnalysisObject::operator =;

    /// @brief Nullary constructor
    Cutflow() : BaseT(), icurr(0) {}

    /// @brief Constructor using binning object
    Cutflow(const BinningT& binning) : BaseT(binning), icurr(0) { }

    /// @brief Constructor using rvalue binning object
    Cutflow(BinningT&& binning) : BaseT(std::move(binning)), icurr(0) { }

    /// Proper constructor using a cutflow name and a list of edges
    ///
    /// @note the constructor will preprend the empty string as a pre-cut step.
    Cutflow(const std::vector<std::string>& edges, const std::string& path = "")
      : BaseT(std::vector<std::string>{""} + edges, path), icurr(0) { }

    /// @brief Fill function with FillType
    ///
    /// @note This method will set the current fill position icurr.
    virtual int fill(FillType&& coords, const double weight = 1.0, const double fraction = 1.0) {
      int pos = BaseT::fill(std::move(coords), weight, fraction);
      icurr = (size_t)pos; // cannot be negative for discrete AOs
      return pos;
    }

    /// @brief Fill function using an explicit coordinate
    //
    /// @note This method will set the current fill position icurr.
    virtual int fill(std::string val, const double weight = 1.0, const double fraction = 1.0) {
      int pos = BaseT::fill({val}, weight, fraction);
      icurr = (size_t)pos; // cannot be negative for discrete AOs
      return pos;
    }

    /// @brief Convenience method to fill the pre-cut step
    virtual int fillinit(const double weight=1.0, const double fraction=1.0) {
      return fill(""s, weight, fraction);
    }

    /// @brief Convenience method to fill the next cut step if @a cutresult is true
    virtual int fillnext(const bool cutresult, const double weight = 1.0, const double fraction = 1.0) {
      if (!cutresult)  return 0;
      const std::string edge = BaseT::bin(++icurr).xEdge();
      return fill(edge, weight, fraction);
    }

    /// @brief Convenience method to fill the next cut step
    virtual int fillnext(const double weight = 1.0, const double fraction = 1.0) {
      return fillnext(true, weight, fraction);
    }

    /// @brief Convenience method to fill the next cut steps from an n-element results vector @a cutresult
    virtual int fillnext(const std::vector<bool>& cutresults, const double weight = 1.0, const double fraction = 1.0) {
      if (icurr+cutresults.size() > BaseT::numBins()+1) {
        throw RangeError("Number of filled cut results needs to match the Cutflow construction (in cutflow '"+BaseT::path()+"')");
      }
      for (size_t i = 0; i < cutresults.size(); ++i) {
        if (!cutresults[i])  return 0;
        fillnext(weight, fraction); // also incremenrs icurr
      }
      return (int)icurr;
    }

    /// Returns current cutflow index for debugging
    size_t currentIndex() const {  return icurr; }

    /// Scale the cutflow weights by the given factor
    void scale(const double factor) {
      BaseT::scaleW(factor);
    }

    /// Scale the cutflow weights so that the weight count after cut @a edge is @a norm
    void normalizeStep(const std::string& edge, const double norm) {
      const auto& b = BaseT::binAt(edge);
      if (b.sumW() == 0) {
        MSG_WARNING("Failed to scale Cutflow " << BaseT::path() << " as bin " << edge << " is empty");
        return;
      }
      scale(norm/b.sumW());
    }

    /// Alias to scale the cutflow weights so that the weight after the pre-init cut is @a norm
    void normalizeFirst(const double norm) {
      normalizeStep(""s, norm);
    }

    /// Create a string representation
    string str() const {
      using namespace std;
      stringstream ss;
      ss << fixed << std::setprecision(1) << BaseT::bin(1).sumW();
      const size_t weight0len = ss.str().length();
      ss << fixed << std::setprecision(1) << BaseT::bin(1).effNumEntries();
      const size_t count0len = ss.str().length();
      ss.str("");
      ss << BaseT::path() << " cut-flow:\n";
      size_t maxnamelen = 0;
      for (const string& edge : BaseT::xEdges()) {
        maxnamelen = std::max(edge.length(), maxnamelen);
      }
      ss << setw(maxnamelen+5) << "" << "   "
         << setw(weight0len) << right << "Weight" << "   "
         << setw(count0len) << right << "Count" << "    "
         << setw(6) << right << "A_cumu" << "    "
         << setw(6) << right << "A_incr";
      //for (size_t i = 0; i <= ncuts; ++i) {
      const double wtot = BaseT::bin(1).sumW();
      double wlast = wtot;
      for (const auto& bin : BaseT::bins()) {
        const size_t i = bin.index();
        const int pcttot = (wtot == 0) ? -1 : round(100*bin.sumW()/wtot);
        const int pctinc = (i == 1 || wlast == 0) ? -1 : round(100*bin.sumW()/wlast);
        wlast = bin.sumW();
        stringstream ss2;
        ss2 << fixed << setprecision(1) << bin.sumW();
        const string weightstr = ss2.str(); ss2.str("");
        ss2 << fixed << setprecision(1) << bin.effNumEntries();
        const string countstr = ss2.str(); ss2.str("");
        ss2 << fixed << setprecision(3) << pcttot << "%";
        const string pcttotstr = ss2.str(); ss2.str("");
        ss2 << fixed << setprecision(3) << pctinc << "%";
        const string pctincstr = ss2.str();
        ss << "\n"
           << setw(maxnamelen+5) << left << (i == 1 ? "" : "Pass "+BaseT::bin(i).xEdge()) << "   "
           << setw(weight0len) << right << weightstr << "    "
           << setw(count0len) << right << countstr << "    "
           << setw(6) << right << (pcttot < 0 ? "- " : pcttotstr) << "    "
           << setw(6) << right << (pctinc < 0 ? "- " : pctincstr);
      }
      return ss.str();
    }

    /// Print string representation to a stream
    void print(std::ostream& os) const {
      os << str() << std::flush;
    }

  protected:

    /// Get a logger object.
    Log& getLog() const {
      return Rivet::Log::getLog("Rivet.Cutflow");
    }

    size_t icurr;

  };


  /// FillCollector specialisation for all Cutflow-like AO
  template <>
  class FillCollector<Cutflow> : public Cutflow {
  public:

    using YAO = Cutflow;
    using Ptr = shared_ptr<FillCollector<YAO>>;
    using YAO::operator =;

    FillCollector() : YAO() { }

    /// Constructor
    ///
    /// We call the cheaper constructor based on the
    /// binning to avoid copying of the bin content.
    /// The underlying binning object is still used
    /// in analize() by many routines, e.g. to query
    /// numBins() or to loop over bins() with
    /// subsequent calls to bin xEdge() etc.
    FillCollector(typename YAO::Ptr yao) : YAO(yao->binning()) {
      YAO::setPath(yao->path());
    }

    /// Overloaded fill method, which stores Fill info
    /// until Multiplexer<T>::collapseEventGroup() is called.
    ///
    /// @todo Do we need to deal with users using fractions directly?
    int fill(typename YAO::FillType&& fillCoords,
             const double weight=1.0, const double fraction=1.0) {
      (void)fraction; // suppress unused variable warning
      YAO::icurr = YAO::_binning.globalIndexAt(fillCoords);
      _fills.insert(_fills.end(), { std::move(fillCoords), weight } );
      return (int)YAO::icurr;
    }

    int fill(const std::string val, const double weight = 1.0, const double fraction = 1.0) {
      return fill(typename YAO::FillType{val}, weight, fraction);
    }

    int fillnext(const bool cutresult, const double weight = 1.0, const double fraction = 1.0) {
      if (!cutresult)  return 0;
      return fill(YAO::bin(++YAO::icurr).xEdge(), weight, fraction);
    }

    int fillnext(const double weight = 1.0, const double fraction = 1.0) {
      return fillnext(true, weight, fraction);
    }

    int fillnext(const std::vector<bool>& cutresults, const double weight = 1.0, const double fraction = 1.0) {
      if (YAO::icurr+cutresults.size() > YAO::numBins()+1) {
        throw RangeError("Number of filled cut results needs to match the Cutflow construction (in cutflow '"+YAO::path()+"')");
      }
      for (size_t i = 0; i < cutresults.size(); ++i) {
        if (!cutresults[i])  return 0;
        fillnext(weight, fraction);
      }
      return (int)YAO::icurr;
    }

    int fillinit(const double weight=1.0, const double fraction=1.0) {
      return fill(""s, weight, fraction);
    }

    /// Empty the subevent stack (for start of new event group).
    void reset() noexcept { _fills.clear(); }

    /// Access the fill info subevent stack.
    const Fills<YAO>& fills() const { return _fills; }

  private:

    Fills<YAO> _fills;

  };

  /// @name Cutflow utilities
  /// @{

  /// Convenience alias
  using CutflowPtr = MultiplexPtr<Multiplexer<Cutflow>>;

  /// Print a Cutflow to a stream
  inline std::ostream& operator << (std::ostream& os, const Cutflow& cf) {
    return os << cf.str();
  }

  inline std::ostream& operator << (std::ostream& os, const CutflowPtr& cf) {
    return os << cf->str();
  }

  /// @}


  /// A container for several Cutflow objects, with some convenient batch access
  class Cutflows : public YODA::FillableStorage<2, CutflowPtr, std::string>,
                     public YODA::Fillable {
  public:

    using BaseT = YODA::FillableStorage<2, CutflowPtr, std::string>;
    using BinContentT = CutflowPtr;
    using BinT = YODA::Bin<1, BinContentT, BaseT::BinningT>;
    using FillDim = std::integral_constant<size_t, 2>;

    /// @name Constructors
    /// @{

    /// Nullary constructor
    Cutflows() : BaseT(groupAdapter<FillDim{}, BinContentT, std::string>) { }

    /// Constructor using a vector of (outer) edges
    Cutflows(const std::vector<std::string>& edges)
         : BaseT(YODA::Axis<std::string>(edges), groupAdapter<FillDim{}, BinContentT, std::string>) { }

    /// Constructor using an initializer_set of (outer) edges
    Cutflows(std::initializer_list<std::string>&& edges)
         : BaseT(YODA::Axis<std::string>(std::move(edges)), groupAdapter<FillDim{}, BinContentT, std::string>) { }

    /// @}

    /// @name Fill methods
    /// @{

    /// Fill method using the FillType of the underlying FIllableStorage
    int fill(const std::string& grpCoord, const std::string& edge, const double weight = 1.0, const double fraction = 1.0) {
      auto& cfl = BaseT::binAt({grpCoord});
      if (!cfl.raw())  return -1; // nullptr if bin not booked
      return cfl->fill({edge}, weight, fraction);
    }

    /// Method to fill the pre-cut steps
    int groupfillinit(const double weight=1.0, const double fraction=1.0) {
      bool pass = true;
      for (auto& cfl : BaseT::bins()) {
        if (!cfl.get())  continue;
        pass &= cfl->fillinit(weight, fraction);
      }
      return (int)pass;
    }

    /// Method to fill the next step if @a cutresult is true
    int groupfillnext(const bool cutresult, const double weight = 1.0, const double fraction = 1.0) {
      if (!cutresult)  return 0;
      bool pass = true;
      for (auto& cfl : BaseT::bins()) {
        if (!cfl.get())  continue;
        pass &= cfl->fillnext(cutresult, weight, fraction);
      }
      return (int)pass;
    }

    /// Short-hand method that assumes the cut result is true
    int groupfillnext(const double weight = 1.0, const double fraction = 1.0) {
      return groupfillnext(true, weight, fraction);
    }

    /// Method using a vector of @a cutresults
    int groupfillnext(const std::vector<bool>& cutresults, const double weight = 1.0, const double fraction = 1.0) {
      bool pass = true;
      for (auto& cfl : BaseT::bins()) {
        if (!cfl.get())  continue;
        pass &= cfl->fillnext(cutresults, weight, fraction);
      }
      return (int)pass;
    }

    /// Method that only fills the next step at coordinate @a grpCoord
    int fillnext(const std::string& grpCoord, const double weight=1.0, const double fraction=1.0) {
      return BaseT::binAt(grpCoord)->fillnext(weight, fraction);
    }

    /// Method that only fills the next step at coordinate @a grpCoord if @a cutresult is true
    int fillnext(const std::string& grpCoord, const bool cutresult, const double weight=1.0, const double fraction=1.0) {
      return BaseT::binAt(grpCoord)->fillnext(cutresult, weight, fraction);
    }

    /// Method that only fills the next step at coordinate @a grpCoord using a vector of @a cutresults
    int fillnext(const std::string& grpCoord, const std::vector<bool>& cutresults,
                 const double weight = 1.0, const double fraction = 1.0) {
      return BaseT::binAt(grpCoord)->fillnext(cutresults, weight, fraction);
    }

    /// @}

    /// @name Reset methods
    /// @{

    /// @brief Reset the histogram.
    ///
    /// Keep the binning but set all bin contents and related quantities to zero
    void reset() noexcept { BaseT::reset(); }

    /// @}

    /// @name Whole histo data
    /// @{

    /// @brief Fill dimension of the object (number of conent axes + temprary axis)
    size_t fillDim() const noexcept { return FillDim{}; }

    /// @brief Total dimension of the object (number of fill axes + filled value)
    size_t dim() const noexcept { return fillDim()+1; }

    /// @brief Get the number of fills (fractional fills are possible).
    double numEntries(const bool includeOverflows=true) const noexcept {
      double n = 0;
      for (const auto& cfl : BaseT::bins(includeOverflows)) {
        if (!cfl.get())  continue;
        n += cfl->numEntries(includeOverflows);
      }
      return n;
    }

    /// @brief Get the effective number of fills.
    double effNumEntries(const bool includeOverflows=true) const noexcept {
      double n = 0;
      for (const auto& cfl : BaseT::bins(includeOverflows)) {
        if (!cfl.get())  continue;
        n += cfl->effNumEntries(includeOverflows);
      }
      return n;
    }

    /// @brief Calculates sum of weights in histo group.
    double sumW(const bool includeOverflows=true) const noexcept {
      double sumw = 0;
      for (const auto& cfl : BaseT::bins(includeOverflows)) {
        if (!cfl.get())  continue;
        sumw += cfl->sumW(includeOverflows);
      }
      return sumw;
    }

    /// @brief Calculates sum of squared weights in histo group.
    double sumW2(const bool includeOverflows=true) const noexcept {
      double sumw2 = 0;
      for (const auto& cfl : BaseT::bins(includeOverflows)) {
        if (!cfl.get())  continue;
        sumw2 += cfl->sumW2(includeOverflows);
      }
      return sumw2;
    }

    /// @brief Get the total volume of the histogram group.
    double integral(const bool includeOverflows=true) const noexcept {
      return sumW(includeOverflows);
    }

    /// @brief Get the total volume error of the histogram group.
    double integralError(const bool includeOverflows=true) const noexcept {
      return sqrt(sumW2(includeOverflows));
    }

    /// @}

    /// @name Transformations
    /// @{

    /// @brief Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(const double scalefactor) noexcept {
      for (auto& cfl : BaseT::bins(true, true)) {
        if (!cfl.get())  continue;
        cfl->scaleW(scalefactor);
      }
    }

    /// @brief Rescale as if all fill weights had been different by factor @a scalefactor along dimension @a i.
    void scale(const double scalefactor) noexcept {
      for (auto& cfl : BaseT::bins(true, true)) {
        if (!cfl.get())  continue;
        cfl->scale(scalefactor);
      }
    }

    /// Scale the cutflow weights so that the weight count after cut @a edge is @a norm
    ///
    /// @note Consider each sub-histogram separately when calculating the integral.
    void normalizeStep(const std::string& edge, const double norm) {
      for (auto& cfl : BaseT::bins(true, true)) {
        if (!cfl.get())  continue;
        if (cfl->binAt(edge).sumW() != 0.)  cfl->normalizeStep(edge, norm);
      }
    }

    /// Alias to scale the cutflow weights so that the weight after the pre-init cut is @a norm
    void normalizeFirst(const double norm) {
      normalizeStep(""s, norm);
    }

    /// @}

    /// @name Streaming utilities
    /// @{

    /// Create a string representation
    string str() const {
      std::stringstream ss;
      for (auto& cf : BaseT::bins(true, true)) {
        if (!cf.get())  continue;
        ss << cf << "\n\n";
      }
      return ss.str();
    }

    /// Print string representation to a stream
    void print(std::ostream& os) const {
      os << str() << std::flush;
    }

    /// @}

  };

  /// @name Cutflows utilities
  /// @{

  /// Convenience alias
  using CutflowsPtr = std::shared_ptr<Cutflows>;

  /// Print a Cutflows to a stream
  inline std::ostream& operator << (std::ostream& os, const Cutflows& cfs) {
    return os << cfs.str();
  }

  /// Print a Cutflows to a stream
  inline std::ostream& operator << (std::ostream& os, const CutflowsPtr& cfs) {
    return os << cfs->str();
  }

  /// @}

}

#endif
