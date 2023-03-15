// -*- C++ -*-
#ifndef RIVET_TYPEREGISTRY_HH
#define RIVET_TYPEREGISTRY_HH

#include "YODA/AnalysisObject.h"
#include "YODA/ReaderYODA.h"

#include <unordered_map>
#include <memory>

namespace YODA {

  using AnalysisObjectPtr = std::shared_ptr<YODA::AnalysisObject>;

}

namespace Rivet {


  /// @brief A polymorphic base type for the AO type handles
  struct TypeBaseHandle {

    TypeBaseHandle() = default;

    virtual ~TypeBaseHandle() { }

    virtual bool copyAO(YODA::AnalysisObjectPtr src,
                        YODA::AnalysisObjectPtr dst,
                        const double scale = 1.0) const = 0;

    virtual bool addAO(YODA::AnalysisObjectPtr src,
                       YODA::AnalysisObjectPtr& dst,
                       const double scale = 1.0) const = 0;

  };



  /// @brief The type-specific handle that can perform
  /// type-specific operations for objects of type T
  template<typename T>
  struct TypeHandle : public TypeBaseHandle {

    bool addAO(YODA::AnalysisObjectPtr src,
               YODA::AnalysisObjectPtr& dst,
               const double scale = 1.0) const;

    bool copyAO(YODA::AnalysisObjectPtr src,
                YODA::AnalysisObjectPtr dst,
                const double scale = 1.0) const;

  };


  /// @brief Simple registry to keep track of AO types known to Rivet
  class TypeRegistry {
  protected:

    /// Some convenient aliases
    using TypeHandlePtr = std::shared_ptr<TypeBaseHandle>;
    using TypeMap = std::unordered_map<string, TypeHandlePtr>;
    using TypeItr = typename TypeMap::const_iterator;

  public:

    /// Singleton creation function
    static TypeRegistry& create();

    /// disable copy and move constructors
    TypeRegistry(const TypeRegistry&) = delete;
    TypeRegistry(TypeRegistry&&) = delete;

    /// disable copy and move assignments
    TypeRegistry& operator=(const TypeRegistry&) = delete;
    TypeRegistry& operator=(TypeRegistry&&) = delete;

    /// Adds type T with AO::type() @a name to the TypeMap
    template<typename T>
    void addType(const string& name);

    /// Returns type with AO::type() @a name from the TypeMap
    TypeHandlePtr findType(const string& name) const;


  private:

    /// Loads a default set of YODA AO types into TypeMap
    template<typename... Args>
    void loadDefaultTypes();

    TypeRegistry() = default; // private since it's a singleton

    /// The TypeMap of all registered types
    TypeMap _register;

    /// A reference to YODA's Reader singleton
    YODA::Reader& _reader = YODA::ReaderYODA::create();

  };


  /// Registers a new type with AO::type() name @a name
  template<typename T>
  void registerType(const string& name);

  /// Retrieves the type handle from the registry
  std::shared_ptr<TypeBaseHandle> findRegisteredType(const string& typeName);


}

#endif
