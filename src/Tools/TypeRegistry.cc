// -*- C++ -*-
#include "Rivet/Tools/TypeRegistry.hh"
#include "Rivet/Tools/TypeTraits.hh"
#include "YODA/Counter.h"
#include "YODA/Histo.h"
#include "YODA/Profile.h"
#include "YODA/Scatter.h"

using namespace std;

namespace {

  /// Fold expression to call lambda @a F with
  /// each argument in the variadic set
  template<typename... Args, typename F>
  constexpr void for_each_arg(F&& f) {
    (( f(Args{}) ), ...);
  }

}

namespace Rivet {

  /// Singleton creation function
  TypeRegistry& TypeRegistry::create() {
    static TypeRegistry _instance;
    _instance.loadDefaultTypes<double,int,string>();
    return _instance;
  }

  /// Registers a new type with AO::type() @a name
  template<typename T>
  void registerType(const string& name) {
    TypeRegistry& reg = TypeRegistry::create();
    reg.addType<T>(name);
  }

  /// Retrieves the type handle from the registry
  std::shared_ptr<TypeBaseHandle> findRegisteredType(const string& typeName) {
    TypeRegistry& reg = TypeRegistry::create();
    return reg.findType(typeName);
  }


  /// If @a src and @a dst both are of same subclass T, add the
  /// contents of @a src into @a dst and return true, otherwise
  /// return false. Apply scaling @a scale in case T is fillable.
  /// If @a dst is nullptr, assign @a src.
  template<typename T>
  bool TypeHandle<T>::addAO(YODA::AnalysisObjectPtr src,
                            YODA::AnalysisObjectPtr& dst,
                            const double scale) const {
    if constexpr (isFillable<T>::value) {
      shared_ptr<T> srcPtr = std::static_pointer_cast<T>(src);
      srcPtr->scaleW(scale);
      if (dst == nullptr) { dst = src; return true; }
      try { *std::static_pointer_cast<T>(dst) += *srcPtr; }
      catch (YODA::BinningError&) { return false; }
      return true;
    }
    else if (dst == nullptr) { dst = src; return true; }
    return false;
  }


  /// If @a src and @a dst both are of same subclass T, copy the
  /// contents of @a src into @a dst and return true, otherwise
  /// return false. Apply scaling @a scale in case T is fillable.
  template<typename T>
  bool TypeHandle<T>::copyAO(YODA::AnalysisObjectPtr src,
                             YODA::AnalysisObjectPtr dst,
                             const double scale) const {
    if (src->type() != dst->type()) {
      throw YODA::LogicError("Operation requries types to be the same!");
    }
    shared_ptr<T> dstPtr = std::static_pointer_cast<T>(dst);
    *dstPtr = *std::static_pointer_cast<T>(src);
    if constexpr (isFillable<T>::value) { dstPtr->scaleW(scale); }
    return true;
  }


  /// If @a src and @a dst both are of same subclass T, add the
  /// contents of @a src to @a dst and return true, otherwise
  /// return false. Apply scaling @a scale in case T is fillable.
  template<typename T>
  void TypeRegistry::addType(const string& name) {
    const TypeItr& res = _register.find(name);
    if (res == _register.end()) {
      _register[name] = make_shared<TypeHandle<T>>();
    }
    _reader.addType<T>(name); // also let YODA know
  }


  /// Check if @a name has been registered and return
  /// associated TypeHandle or nullptr otherwise.
  shared_ptr<TypeBaseHandle> TypeRegistry::findType(const string& name) const {
    const TypeItr& tmp = _register.find(name);
    if (tmp == _register.end())  return nullptr;
    return tmp->second;
  }


  /// Register a default set of types in the calling singleton
  ///
  /// @todo Lots of overlap with a similar piece of logic in YODA
  /// here - can this perhaps be merged...?
  template<typename... Args>
  void TypeRegistry::loadDefaultTypes() {

    // load 0D types
    addType<YODA::Counter>("Counter");
    addType<YODA::Estimate0D>("Estimate0D");

    // load scatters and short-hand types
    auto addShorthands = [&](auto I) {
      addType<YODA::HistoND<I+1>>("Histo"+to_string(I+1)+"D");
      addType<YODA::ProfileND<I+1>>("Profile"+to_string(I+1)+"D");
      addType<YODA::ScatterND<I+1>>("Scatter"+to_string(I+1)+"D");
    };
    MetaUtils::staticFor<3>(addShorthands); // exec for up to 3 dimensions

    // load BinnedHisto/BinnedProfile in 1D
    for_each_arg<Args...>([&](auto&& arg) {
      using A1 = decay_t<decltype(arg)>;
      using BH = YODA::BinnedHisto<A1>;
      addType<BH>(BH().type());
      using BP = YODA::BinnedProfile<A1>;
      addType<BP>(BP().type());
      using BE = YODA::BinnedEstimate<A1>;
      addType<BE>(BE().type());

      // load BinnedHisto/BinnedProfile in 2D
      for_each_arg<Args...>([&](auto&& arg) {
        using A2 = decay_t<decltype(arg)>;
        using BH = YODA::BinnedHisto<A1,A2>;
        addType<BH>(BH().type());
        using BP = YODA::BinnedProfile<A1,A2>;
        addType<BP>(BP().type());
        using BE = YODA::BinnedEstimate<A1,A2>;
        addType<BE>(BE().type());

        // load BinnedHisto/BinnedProfile in 3D
        for_each_arg<Args...>([&](auto&& arg) {
          using A3 = decay_t<decltype(arg)>;
          using BH = YODA::BinnedHisto<A1,A2,A3>;
          addType<BH>(BH().type());
          //using BP = YODA::BinnedProfile<A1,A2,A3>;
          //addType<BP>(BP().type());
          using BE = YODA::BinnedEstimate<A1,A2,A3>;
          addType<BE>(BE().type());
        });
      });
    });

  }

}
