/* boost random/faure.hpp header file
 *
 * Copyright Justinas Vygintas Daugmaudis 2010
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_RANDOM_FAURE_HPP
#define BOOST_RANDOM_FAURE_HPP

#include <boost/random/detail/quasi_random_number_generator_base.hpp>

#include <vector>
#include <cmath>

#include <boost/integer.hpp>

#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <boost/mpl/vector/vector50_c.hpp>
#include <boost/mpl/lower_bound.hpp>

//!\file
//!Describes the quasi-random number generator class template faure.

namespace boost {
namespace random {

/** @cond */
namespace detail {
namespace fr {

struct prime_list
{
  typedef mpl::vector50_c<int,
      2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
     31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
     73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
    127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
    179,  181,  191,  193,  197,  199,  211,  223,  227,  229> type;
};


// Returns the smallest prime greater than or equal to N.
template<int N>
struct prime_ge
{
  typedef prime_list::type prime_list_t;
  typedef typename mpl::lower_bound<prime_list_t, mpl::int_<N> >::type iter_t;

  BOOST_STATIC_CONSTANT(int, value = mpl::deref<iter_t>::type::value);
};

// Returns the integer part of the logarithm base Base of arg.
// In erroneous situations, e.g., integer_log<Base>(0) the function
// returns 0 and does not report the error. This is the intended
// behavior.
template<std::size_t Base>
inline std::size_t integer_log(std::size_t arg)
{
  BOOST_STATIC_ASSERT( 2 <= Base );

  std::size_t ilog = 0;
  while( Base <= arg )
  {
    arg /= Base; ++ilog;
  }
  return ilog;
}

inline std::size_t integer_pow(std::size_t arg, std::size_t p)
{
  std::size_t ipow = 1;
  while( p-- )
    ipow *= arg;
  return ipow;
}

// Computes a table of binomial coefficients modulo qs.
template<typename RealType, std::size_t Dimension>
struct binomial_coefficients
{
  typedef RealType result_type;

  BOOST_STATIC_CONSTANT(std::size_t, dimension_value = Dimension);
  BOOST_STATIC_CONSTANT(int, qs_base = prime_ge<Dimension>::value);

  // The maximum prime value that fits into vector50_c is 229, which
  // means that our primes fit into the uint8_t range.
  // This merits our attention, because binomial values modulo qs_base
  // will never be bigger than qs_base. We can choose an appropriate
  // integer type to hold modulo values and shave off memory footprint.

  typedef typename uint_value_t<qs_base>::fast packed_uint_t;

  inline static RealType inv_qs_base()
  {
    return static_cast<RealType>(1) / static_cast<RealType>(qs_base);
  }

  // default copy c-tor is fine

  binomial_coefficients(std::size_t seq) // c-tor to initialize the binomial table
  {
    resize(n_elements(seq));
  }

  void update(std::size_t seq, RealType (&quasi)[Dimension])
  {
    const std::size_t hisum = n_elements(seq);
    if( coeff.size() != size_hint(hisum) )
      resize(hisum);

    quasi[0] =
      compute_recip(seq, hisum, ytemp.rbegin(), ytemp.rend());

    // Find components using the Faure method.
    for( std::size_t k = 1; k < Dimension; ++k )
    {
      quasi[k] = RealType();
      RealType r = inv_qs_base();

      for( std::size_t i = 0; i != hisum; ++i )
      {
        RealType ztemp = RealType();
        for( std::size_t j = i; j != hisum; ++j )
          ztemp += ytemp[j] * upper_element(i, j, hisum);

        // Sum ( J <= I <= HISUM ) ( old ytemp(i) * binom(i,j) ) mod QS.
        ytemp[i] = std::fmod(ztemp, static_cast<RealType>(qs_base));
        quasi[k] += ytemp[i] * r;
        r *= inv_qs_base();
      }
    }
  }

  void incremental_update(std::size_t& v, RealType (&quasi)[Dimension])
  {
    std::size_t old_v = v++;
    if( old_v >= v ) // overflow in v
      boost::throw_exception( std::overflow_error("faure::checked_increment") );
    else
      update(v, quasi);
  }

private:
  inline static std::size_t n_elements(std::size_t seq)
  {
    return integer_log<qs_base>(seq) + 1;
  }

  inline static std::size_t size_hint(std::size_t n)
  {
    return n * (n + 1) / 2;
  }

  packed_uint_t& upper_element(std::size_t i, std::size_t j, std::size_t dim)
  {
    BOOST_ASSERT( i < dim );
    BOOST_ASSERT( j < dim );
    BOOST_ASSERT( i <= j );
    return coeff[(i * (2 * dim - i + 1)) / 2 + j - i];
  }

  template<typename Iterator>
  inline static RealType compute_recip(std::size_t seq, std::size_t n,
      Iterator first, Iterator last)
  {
    //BOOST_ASSERT( std::distance(first, last) == n );

    // Here we do
    //   Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS**J
    //   Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
    // in one go
    RealType v, r = RealType();
    std::size_t m, k = integer_pow(qs_base, n - 1);
    for( ; first != last; ++first, seq = m, k /= qs_base )
    {
      m  = seq % k;
      v  = (seq - m) / k;
      r += v;
      r *= inv_qs_base();
      *first = v; // saves double dereference
    }
    return r;
  }

  void compute_coefficients(const std::size_t n)
  {
    // Resize and initialize to zero
    coeff.resize(size_hint(n));
    std::fill(coeff.begin(), coeff.end(), 0);

    // The first row and the diagonal is assigned to 1
    upper_element(0, 0, n) = 1;
    for( std::size_t i = 1; i < n; ++i )
    {
      upper_element(0, i, n) = 1;
      upper_element(i, i, n) = 1;
    }

    // Computes binomial coefficients MOD qs_base
    for( std::size_t i = 1; i < n; ++i )
    {
      for( std::size_t j = i + 1; j < n; ++j )
      {
        upper_element(i, j, n) = ( upper_element(i, j-1, n) +
                                   upper_element(i-1, j-1, n) ) % qs_base;
      }
    }
  }

  void resize(std::size_t n)
  {
    ytemp.resize(n);
    compute_coefficients(n);
  }

private:
  // here we cache precomputed data; note that binomial coefficients have
  // to be recomputed iff the integer part of the logarithm of seq changes,
  // which happens relatively rarely.
  std::vector<packed_uint_t> coeff; // packed upper (!) triangular matrix
  std::vector<RealType> ytemp;
};

}} // namespace detail::fr
/** @endcond */

//!class template faure implements a quasi-random number generator as described in
//! \blockquote
//!Henri Faure,
//!Discrepance de suites associees a un systeme de numeration (en dimension s),
//!Acta Arithmetica,
//!Volume 41, 1982, pages 337-351.
//! \endblockquote
//
//! \blockquote
//!Bennett Fox,
//!Algorithm 647:
//!Implementation and Relative Efficiency of Quasirandom
//!Sequence Generators,
//!ACM Transactions on Mathematical Software,
//!Volume 12, Number 4, December 1986, pages 362-376.
//! \endblockquote
//!
//!\attention\b Important: This implementation supports up to 229 dimensions.
//!
//!In the following documentation @c X denotes the concrete class of the template
//!faure returning objects of type @c RealType, u and v are the values of @c X.
//!
//!Some member functions may throw exceptions of type @c std::bad_alloc.
//!
//! \copydoc friendfunctions
template<typename RealType, std::size_t Dimension>
class faure : public detail::quasi_random_number_generator_base<
                        faure<RealType, Dimension>
                      , detail::fr::binomial_coefficients<RealType, Dimension>
                      >
{
/** @cond */
#if defined(BOOST_NO_STATIC_ASSERT)
  BOOST_STATIC_ASSERT( Dimension <= 229 );
#else
  static_assert(Dimension <= 229, "The Faure quasi-random number generator only supports 229 dimensions.");
#endif
/** @endcond */

  typedef faure<RealType, Dimension> self_t;
  typedef detail::fr::binomial_coefficients<RealType, Dimension> lattice_t;
  typedef detail::quasi_random_number_generator_base<self_t, lattice_t> base_t;

  friend class detail::quasi_random_number_generator_base<
                          self_t
                        , lattice_t >;
public:
  typedef RealType result_type;

  /** @copydoc boost::random::niederreiter_base2::min() */
  static result_type min /** @cond */ BOOST_PREVENT_MACRO_SUBSTITUTION /** @endcond */ () { return static_cast<RealType>(0); }

  /** @copydoc boost::random::niederreiter_base2::max() */
  static result_type max /** @cond */ BOOST_PREVENT_MACRO_SUBSTITUTION /** @endcond */ () { return static_cast<RealType>(1); }

  /** @copydoc boost::random::niederreiter_base2::dimension() */
  static std::size_t dimension() { return Dimension; }

  //!Effects: Constructs the default Faure quasi-random number generator.
  //!
  //!Throws: bad_alloc.
  faure()
    : base_t(0) // initialize the binomial table here
  {}

  //!Effects: Constructs the Faure quasi-random number generator,
  //!equivalent to X u; u.seed(init).
  //!
  //!Throws: bad_alloc.
  explicit faure(std::size_t init)
    : base_t(init) // initialize the binomial table here
  {}

  /** @copydetails boost::random::niederreiter_base2::seed()
   * Throws: bad_alloc.
   */
  void seed()
  {
    seed(0);
  }

  /** @copydetails boost::random::niederreiter_base2::seed(std::size_t)
   * Throws: bad_alloc.
   */
  void seed(std::size_t init)
  {
    this->curr_elem = 0;
    this->seq_count = init;
    this->lattice.update(this->seq_count, this->quasi_state);
  }

  //=========================Doxygen needs this!==============================

  //!Requirements: *this is mutable.
  //!
  //!Returns: Returns a successive element of an s-dimensional
  //!(s = X::dimension()) vector at each invocation. When all elements are
  //!exhausted, X::operator() begins anew with the starting element of a
  //!subsequent s-dimensional vector.
  //!
  //!Throws: bad_alloc.

  // Fixed in Doxygen 1.7.0 -- id 612458: Fixed problem handling @copydoc for function operators.
  result_type operator()()
  {
    return base_t::operator()();
  }

  /** @copydoc boost::random::niederreiter_base2::discard(std::size_t)
   * Throws: bad_alloc.
   */
  void discard(std::size_t z)
  {
    base_t::discard(z);
  }

private:
/** @cond hide_private_members */
  void compute_next_vector()
  {
    this->lattice.incremental_update(this->seq_count, this->quasi_state);
  }
/** @endcond */
};

} // namespace random

//It would have been nice to do something like this, but it seems that doxygen
//simply won't show those typedefs, so we spell them out by hand.

/*
#define BOOST_FAURE_TYPEDEF(z, N, text) \
typedef random::faure<double, N> faure_##N##d; \
//

BOOST_PP_REPEAT_FROM_TO(1, 229, BOOST_FAURE_TYPEDEF, unused)

#undef BOOST_FAURE_TYPEDEF
*/

typedef random::faure<double, 1> faure_1d;
typedef random::faure<double, 2> faure_2d;
typedef random::faure<double, 3> faure_3d;
typedef random::faure<double, 4> faure_4d;
typedef random::faure<double, 5> faure_5d;
typedef random::faure<double, 6> faure_6d;
typedef random::faure<double, 7> faure_7d;
typedef random::faure<double, 8> faure_8d;
typedef random::faure<double, 9> faure_9d;
typedef random::faure<double, 10> faure_10d;
typedef random::faure<double, 11> faure_11d;
typedef random::faure<double, 12> faure_12d;

} // namespace boost

#endif // BOOST_RANDOM_FAURE_HPP
