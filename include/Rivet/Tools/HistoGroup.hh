// -*- C++ -*-
#ifndef RIVET_HISTOGROUP_HH
#define RIVET_HISTOGROUP_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {

  namespace {

    /// @brief Dummy fill adapter for the HistoGroup
    ///
    /// The group bin is already being determined in the fill method.
    template <size_t FillDim, typename T, typename AxisT>
    typename YODA::FillableStorage<FillDim, T, AxisT>::FillAdapterT
    groupAdapter = [](auto& /* ptr */, auto&& /* coords */, double /* weight */, double /* fraction */) { };

  }

  template <typename GroupAxisT, typename... AxisT>
  class HistoGroup : public YODA::FillableStorage<sizeof...(AxisT)+1, BinnedHistoPtr<AxisT...>, GroupAxisT>,
                     public YODA::Fillable {
  public:

    using BaseT = YODA::FillableStorage<sizeof...(AxisT)+1, BinnedHistoPtr<AxisT...>, GroupAxisT>;
    using BinContentT = BinnedHistoPtr<AxisT...>;
    using FillDim = std::integral_constant<size_t, sizeof...(AxisT)+1>;

    HistoGroup() : BaseT(groupAdapter<FillDim{}, BinContentT, GroupAxisT>) { }

    HistoGroup(const std::vector<GroupAxisT>& edges)
         : BaseT(YODA::Axis<GroupAxisT>(edges), groupAdapter<FillDim{}, BinContentT, GroupAxisT>) { }

    HistoGroup(std::initializer_list<GroupAxisT>&& edges)
         : BaseT(YODA::Axis<GroupAxisT>(std::move(edges)), groupAdapter<FillDim{}, BinContentT, GroupAxisT>) { }

    int fill(GroupAxisT tmpCoord, AxisT... coords, const double weight = 1.0, const double fraction = 1.0) {
      auto& bin = BaseT::binAt(tmpCoord);
      if (!bin.raw())  return -1; // nullptr if bin not booked
      return bin->fill({coords...}, weight, fraction);
    }

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
      for (const auto& b : BaseT::bins(includeOverflows)) {
        if (!b.get())  continue;
        n += b->numEntries(includeOverflows);
      }
      return n;
    }

    /// @brief Get the effective number of fills.
    double effNumEntries(const bool includeOverflows=true) const noexcept {
      double n = 0;
      for (const auto& b : BaseT::bins(includeOverflows)) {
        if (!b.get())  continue;
        n += b->effNumEntries(includeOverflows);
      }
      return n;
    }

    /// @brief Calculates sum of weights in histo group.
    double sumW(const bool includeOverflows=true) const noexcept {
      double sumw = 0;
      for (const auto& b : BaseT::bins(includeOverflows)) {
        if (!b.get())  continue;
        sumw += b->sumW(includeOverflows);
      }
      return sumw;
    }

    /// @brief Calculates sum of squared weights in histo group.
    double sumW2(const bool includeOverflows=true) const noexcept {
      double sumw2 = 0;
      for (const auto& b : BaseT::bins(includeOverflows)) {
        if (!b.get())  continue;
        sumw2 += b->sumW2(includeOverflows);
      }
      return sumw2;
    }

    /// @brief Return the vector of sum of weights for each histo in the group.
    vector<double> sumWGroup(const bool includeOverflows=true) const noexcept {
      vector<double> rtn;
      rtn.reserve(BaseT::numBins(true));
      for (const auto& b : BaseT::bins(true)) {
        if (!b.get()) {
          rtn.push_back(0.0);
          continue;
        }
        rtn.push_back( b->sumW(includeOverflows) );
      }
      return rtn;
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

    void divByGroupWidth() const noexcept {
      for (auto& bin : BaseT::bins(true, true)) {
        if (!bin.get())  continue;
        const double bw = bin.dVol();
        if (bw)  bin->scaleW(1.0/bw);
      }
    }

    /// @brief Rescale as if all fill weights had been different by factor @a scalefactor.
    void scaleW(const double scalefactor) noexcept {
      for (auto& bin : BaseT::bins(true, true)) {
        if (!bin.get())  continue;
        bin->scaleW(scalefactor);
      }
    }

    /// @brief Rescale as if all fill weights had been different by factor @a scalefactor along dimension @a i.
    void scale(const size_t i, const double scalefactor) noexcept {
      for (auto& bin : BaseT::bins(true, true)) {
        if (!bin.get())  continue;
        bin->scale(i, scalefactor);
      }
    }

    /// @brief Normalize the (visible) histo "volume" to the @a normto value.
    ///
    /// @note Consider each sub-histogram separately when calculating the integral.
    void normalize(const double normto=1.0, const bool includeOverflows=true) {
      for (auto& bin : BaseT::bins(true, true)) {
        if (!bin.get())  continue;
        if (bin->integral(includeOverflows) != 0.)  bin->normalize(normto, includeOverflows);
      }
    }

    /// @brief Normalize the (visible) histo "volume" to the @a normto value.
    ///
    /// @note Treat the histogram group as an entity when calculating the integral.
    void normalizeGroup(const double normto=1.0, const bool includeOverflows=true) {
      const double oldintegral = integral(includeOverflows);
      if (oldintegral == 0) throw WeightError("Attempted to normalize a histogram group with null area");
      scaleW(normto / oldintegral);
    }

    /// @}

  };

  /// @name Convenient aliases and dimension-specific short-hands
  /// @{

  template <typename GroupAxisT, typename... AxisT>
  using HistoGroupPtr = std::shared_ptr<HistoGroup<GroupAxisT, AxisT...>>;
  using Histo1DGroupPtr = HistoGroupPtr<double,double>;
  using Histo2DGroupPtr = HistoGroupPtr<double,double,double>;

  /// @}

}

#endif
