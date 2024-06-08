// -*- C++ -*-
#ifndef RIVET_TypeTraits_HH
#define RIVET_TypeTraits_HH

#include <type_traits>
#include <iterator>

namespace Rivet {


  /// Mechanisms to allow references and pointers to templated types
  /// to be distinguished from one another (since C++ doesn't allow
  /// partial template specialisation for functions.
  ///
  /// Traits methods use specialisation of class/struct templates, and
  /// some trickery with typedefs and static const integral types (or
  /// enums) to implement partial function specialisation as a work-around.

  /// @cond INTERNAL

  struct RefType { };

  struct PtrType { };

  template <typename T>
  struct TypeTraits;

  template <typename U>
  struct TypeTraits<const U&> {
    using ArgType = RefType;
  };

  template <typename U>
  struct TypeTraits<const U*> {
    using ArgType = PtrType;
  };

  template<typename T>
  struct TypeIdentity { using type = T; };


  /// Conditional to return the signed integer type, otherwise same type
  template<typename T>
  using sign_if_integral = std::conditional_t<(std::is_integral_v<T> && !std::is_same_v<T,bool>),
                                        std::make_signed<T>, TypeIdentity<T>>;
  //
  template<typename T>
  using sign_if_integral_t = typename sign_if_integral<T>::type;

  /// Conditional to return true if all input parameters are unsgined
  template<typename... Args>
  using all_unsigned = std::conjunction<std::is_unsigned<Args>...>;

  /// Conditional to return signed common type if input types are mixed signed/unsigned
  template<typename... Args>
  using signed_if_mixed_t = std::conditional_t<all_unsigned<Args...>::value,
                                               std::common_type_t<Args...>,
                                               std::common_type_t<sign_if_integral_t<Args>...>>;


  /// SFINAE definition of dereferenceability trait, cf. Boost has_dereference
  template <typename T, typename=void>
  struct Derefable : std::false_type {};
  //
  template <typename T>
  struct Derefable<T, std::void_t< decltype(*std::declval<T>())> > : std::true_type {};


  /// SFINAE struct to check for iterator concept
  ///
  /// @note Also works for C-style arrays.
  template<typename T, typename = void>
  struct Iterable : std::false_type { };
  //
  template<typename T>
  struct Iterable<T, std::void_t<std::decay_t<decltype(std::begin(std::declval<const T&>()))>,
                                 std::decay_t<decltype(std::end(std::declval<const T&>()))>>>  : std::true_type { };

  template<typename T>
  inline constexpr bool is_iterable_v = Iterable<T>::value;

  //template<typename... T>
  //inline constexpr bool isIterable = std::conjunction<Iterable<T>...>::value;
  template<typename T>
  using isIterable = std::enable_if_t<Iterable<T>::value>;


  // SFINAE struct to check for const_iterator concept
  template<typename T, typename = void>
  struct CIterable : std::false_type { };

  template<typename T>
  struct CIterable<T, std::void_t<decltype(std::declval<typename std::decay_t<T>::const_iterator>())>> : std::true_type { };

  template<typename T>
  inline constexpr bool is_citerable_v = CIterable<T>::value;

  //template<typename... T>
  //inline constexpr bool isCIterable = std::conjunction<CIterable<T>...>::value;
  template<typename... Args>
  using isCIterable = std::enable_if_t<(CIterable<Args>::value && ...)>;

  // SFINAE struct to check for cstr
  template<typename T>
  struct isCString : std::false_type { };

  template<>
  struct isCString<char[]> : std::true_type { };

  template<size_t N>
  struct isCString<char[N]> : std::true_type { };

  template <typename T>
  inline constexpr bool is_cstring_v = isCString<T>::value;


  /// SFINAE check if T has a binning() method
  template<typename T, typename = void>
  struct hasBinning : std::false_type { };
  //
  template<typename T>
  struct hasBinning<T, std::void_t<decltype(std::declval<T>().binning())>> : std::true_type { };


  /// SFINAE check if T is a YODA Fillable
  template<typename T, typename = void>
  struct isFillable : std::false_type { };
  //
  template<typename T>
  struct isFillable<T, std::void_t<decltype(typename T::FillType{})>> : std::true_type { };


  // SFINAE struct to check if U can be used as an argument of T
  template<typename T, typename U, typename = void>
  struct isArgumentOf : std::false_type { };
  //
  template<typename T, typename U>
  struct isArgumentOf<T, U, std::void_t<decltype(T(std::declval<U>()))>>
                  : std::true_type { };

  /// Check if all elements in parameter pack Us can be used as an argument of T
  template <typename T, typename... Us>
  using allArgumentsOf = typename std::conjunction<isArgumentOf<T, Us>...>;

  /// SFINAE check if T has XYZ methods
  template <typename T, typename=void>
  struct HasXYZ : std::false_type {};
  template <typename T>
  struct HasXYZ<T, std::void_t< decltype(std::declval<T>().x() + std::declval<T>().y() +
                                         std::declval<T>().z())> > : std::true_type {};


  /// SFINAE check if T has XYZT methods
  template <typename T, typename=void>
  struct HasXYZT : std::false_type {};
  template <typename T>
  struct HasXYZT<T, std::void_t< decltype(std::declval<T>().x() +  std::declval<T>().y() +
                                          std::declval<T>().z() +  std::declval<T>().t())> > : std::true_type {};


  /// @endcond

}

#endif
