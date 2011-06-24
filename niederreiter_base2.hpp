/* boost random/nierderreiter_base2.hpp header file
 *
 * Copyright Justinas Vygintas Daugmaudis 2010
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_RANDOM_NIEDERREITER_BASE2_HPP
#define BOOST_RANDOM_NIEDERREITER_BASE2_HPP

#include <boost/random/detail/gray_coded_qrng_base.hpp>
#include <boost/random/detail/config.hpp>
#include <boost/random/detail/operators.hpp>

#include <limits>
#include <bitset>
#include <boost/integer/static_log2.hpp>

#include <boost/static_assert.hpp>

#include <boost/cstdint.hpp>

//!\file
//!Describes the quasi-random number generator class template niederreiter_base2.
//!
//!\b Note: it is especially useful in conjunction with class template uniform_real.

namespace boost {
namespace random {

/** @cond */
namespace detail {
namespace nb2 {

/*
  Primitive polynomials in binary encoding
  {
    { 1, 0, 0, 0, 0, 0 },    1
    { 0, 1, 0, 0, 0, 0 },    x
    { 1, 1, 0, 0, 0, 0 },    1 + x
    { 1, 1, 1, 0, 0, 0 },    1 + x + x^2
    { 1, 1, 0, 1, 0, 0 },    1 + x + x^3
    { 1, 0, 1, 1, 0, 0 },    1 + x^2 + x^3
    { 1, 1, 0, 0, 1, 0 },    1 + x + x^4
    { 1, 0, 0, 1, 1, 0 },    1 + x^3 + x^4
    { 1, 1, 1, 1, 1, 0 },    1 + x + x^2 + x^3 + x^4
    { 1, 0, 1, 0, 0, 1 },    1 + x^2 + x^5
    { 1, 0, 0, 1, 0, 1 },    1 + x^3 + x^5
    { 1, 1, 1, 1, 0, 1 },    1 + x + x^2 + x^3 + x^5
    { 1, 1, 1, 0, 1, 1 }     1 + x + x^2 + x^4 + x^5
  };
*/

template<std::size_t D>
struct primitive_polynomial;

#define BOOST_NIEDERREITER_BASE2_PRIM_POLY(D, V) \
  template<> \
   struct primitive_polynomial<D> { \
     static const int value = V; \
     static const int degree = ::boost::static_log2<V>::value; \
   } \
/**/

BOOST_NIEDERREITER_BASE2_PRIM_POLY(0, 1);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(1, 2);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(2, 3);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(3, 7);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(4, 11);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(5, 13);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(6, 19);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(7, 25);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(8, 31);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(9, 37);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(10, 41);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(11, 47);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(12, 55);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(13, 59);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(14, 61);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(15, 67);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(16, 73);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(17, 87);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(18, 91);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(19, 97);
BOOST_NIEDERREITER_BASE2_PRIM_POLY(20, 103);

#undef BOOST_NIEDERREITER_BASE2_PRIM_POLY

// Return the base 2 logarithm for a given bitset v
template<typename Bitset>
inline int ilog2_bitset_trivial(Bitset v)
{
  int ret = 0;
  for( ; v.any(); v >>= 1 )
    ++ret;
  return ret - 1;
}

template<typename Bitset>
inline int ilog2_bitset(Bitset v)
{
  Bitset a = v >> (v.size() / 2);
  return a.any() ? ilog2_bitset(a) + (v.size() / 2):
                   ilog2_bitset_trivial(v);
}


// Multiply polynomials over Z_2.
template<typename Bitset>
inline Bitset modulo2_multiply(int P, Bitset v)
{
  int pos = 0;
  Bitset pt = 0;
  for( ; P; P >>= 1, ++pos )
    if( P & 1 ) pt ^= v << pos;
  return pt;
}


// Calculate the values of the constants V(J,R) as
// described in BFN section 3.3.
//
// px = appropriate irreducible polynomial for current dimension
// pb = polynomial defined in section 2.3 of BFN.
// pb is modified
template<typename Bitset, typename T, int MaxE>
inline void calculate_v(int px, Bitset& pb, int& pb_degree, T (&v)[MaxE])
{
  const T arbitrary_element = 1;  // arbitray element of Z_2

  // Now choose a value of Kj as defined in section 3.3.
  // We must have 0 <= Kj < E*J = M.
  // The limit condition on Kj does not seem very relevant
  // in this program.
  int kj = pb_degree;

  // Now multiply B by PX so B becomes PX**J.
  // In section 2.3, the values of Bi are defined with a minus sign :
  // don't forget this if you use them later!
  pb = modulo2_multiply(px, pb);
  pb_degree = ilog2_bitset(pb);

  // Now choose values of V in accordance with
  // the conditions in section 3.3.
  std::fill(v, v + kj, 0);

  // Quoting from BFN: "Our program currently sets each K_q
  // equal to eq. This has the effect of setting all unrestricted
  // values of v to 1."
  // Actually, it sets them to the arbitrary chosen value.
  // Whatever.
  for( int r = kj; r < pb_degree; ++r) {
    v[r] = arbitrary_element;
  }

  // Calculate the remaining V's using the recursion of section 2.3,
  // remembering that the B's have the opposite sign.
  for(int r = 0; r < MaxE - pb_degree; ++r) {
    T term = 0;
    Bitset pb_c = pb;
    for(int k = 0; k < pb_degree; ++k, pb_c >>= 1)
    {
      if( pb_c.test(0) )
        term ^= v[r + k];
    }
    v[r+pb_degree] = term;
  }
}


template<int Degree, int BitCount>
struct product_degree
{
  BOOST_STATIC_CONSTANT(int, value = Degree * ( (BitCount / Degree) + 1 ));
};

template<std::size_t D>
struct compute_lattice
{
  template<typename T, int BitCount, int MaxE, std::size_t Dimension>
  void operator()(T (&ci)[BitCount][BitCount], T (&v)[MaxE],
                  T (&cj)[BitCount][Dimension]) const
  {
    enum {
      iteration  = Dimension - D,
      poly_index = iteration + 1,
      px_value   = primitive_polynomial<poly_index>::value,
      px_degree  = primitive_polynomial<poly_index>::degree,
      max_degree = product_degree<px_degree, BitCount>::value
    };

    // For each dimension, we need to calculate powers of an
    // appropriate irreducible polynomial, see Niederreiter
    // page 65, just below equation (19).
    // Copy the appropriate irreducible polynomial into PX,
    // and its degree into E.  Set polynomial B = PX ** 0 = 1.
    // M is the degree of B.  Subsequently B will hold higher
    // powers of PX.
    int pb_degree = 0;
    std::bitset<max_degree> pb_value = 1;

    int j = 0;
    while( j < BitCount )
    {
      // If U = 0, we need to set B to the next power of PX
      // and recalculate V.
      calculate_v(px_value, pb_value, pb_degree, v);

      // Niederreiter (page 56, after equation (7), defines two
      // variables Q and U.  We do not need Q explicitly, but we
      // do need U.

      // Advance Niederreiter's state variables.
      for( int u = 0; u < px_degree && j < BitCount; ++u, ++j )
      {
        // Now C is obtained from V.  Niederreiter
        // obtains A from V (page 65, near the bottom), and then gets
        // C from A (page 56, equation (7)).  However this can be done
        // in one step.  Here CI(J,R) corresponds to
        // Niederreiter's C(I,J,R).
        for(int r = 0; r < BitCount; ++r) {
          ci[r][j] = v[r + u];
        }
      }
    }

    // The array CI now holds the values of C(I,J,R) for this value
    // of I.  We pack them into array CJ so that CJ(I,R) holds all
    // the values of C(I,J,R) for J from 1 to NBITS.
    for(int r = 0; r < BitCount; ++r)
    {
      T term = 0;
      for(int j = 0; j < BitCount; ++j)
        term = 2*term + ci[r][j];
      cj[r][iteration] = term;
    }

    compute_lattice<D - 1> ncj; ncj(ci, v, cj);
  }
};

template<>
struct compute_lattice<0>
{
  template<typename T, int BitCount, int MaxE, std::size_t Dimension>
  void operator()(T (&)[BitCount][BitCount], T (&)[MaxE],
                  T (&)[BitCount][Dimension]) const
  {
    // recursion stop
  }
};

} // namespace nb2

template<typename IntType, std::size_t Dimension>
struct niederreiter_base2_lattice
{
#if defined(BOOST_NO_STATIC_ASSERT)
  BOOST_STATIC_ASSERT( Dimension <= 20 );
#else
  static_assert(Dimension <= 20, "The Niederreiter base 2 quasi-random number generator only supports 20 dimensions.");
#endif

  typedef IntType result_type;

  BOOST_STATIC_CONSTANT(std::size_t, dimension_value = Dimension);

  BOOST_STATIC_CONSTANT(int, bit_count = std::numeric_limits<IntType>::digits);
  BOOST_STATIC_CONSTANT(int, prim_degree = nb2::primitive_polynomial<Dimension>::degree);

  // Max degree of a polynomial product that will be computed by modulo2_multiply
  BOOST_STATIC_CONSTANT(int, max_degree = (nb2::product_degree<prim_degree, bit_count>::value));

  // default copy c-tor is fine

  niederreiter_base2_lattice(std::size_t) // c-tor to initialize the lattice
  {
    // initial lattice computation
    IntType ci[bit_count][bit_count];
    IntType v[prim_degree + max_degree];
    nb2::compute_lattice<Dimension> compute; compute(ci, v, bits);
  }

  result_type operator()(int i, int j) const
  {
    return bits[i][j];
  }

  // returns post-increment value
  static std::size_t checked_increment(std::size_t& v)
  {
    std::size_t old_v = v++;
    if( old_v >= v ) // overflow in v
      boost::throw_exception( std::overflow_error("niederreiter_base2::checked_increment") );
    return old_v;
  }

private:
  IntType bits[bit_count][Dimension];
};

} // namespace detail
/** @endcond */

//!class template niederreiter_base2 implements a quasi-random number generator as described in
//! \blockquote
//!Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992).
//! \endblockquote
//!
//!\attention \b Important: This implementation supports up to 20 dimensions.
//!
//!In the following documentation @c X denotes the concrete class of the template
//!niederreiter_base2 returning objects of type @c IntType, u and v are the values of @c X.
//!
//!Some member functions may throw exceptions of type std::overflow_error. This
//!happens when the quasi-random domain is exhausted and the generator cannot produce
//!any more values. The length of the low discrepancy sequence is given by
//! \f$L=Dimension \times 2^{digits}\f$, where digits = std::numeric_limits<IntType>::digits.
template<typename IntType, std::size_t Dimension, IntType c, IntType m>
class niederreiter_base2 : public detail::gray_coded_qrng_base<
                                    niederreiter_base2<IntType, Dimension, c, m>,
                                    detail::niederreiter_base2_lattice<IntType, Dimension> >
{
  typedef niederreiter_base2<IntType, Dimension, c, m> self_t;
  typedef detail::niederreiter_base2_lattice<IntType, Dimension> lattice_t;
  typedef detail::gray_coded_qrng_base<self_t, lattice_t> base_t;

public:
  typedef IntType result_type;

  //!Returns: Tight lower bound on the set of values returned by operator().
  //!
  //!Throws: nothing.
  static result_type min BOOST_PREVENT_MACRO_SUBSTITUTION () { return c == 0u ? 1u: 0u; }

  //!Returns: Tight upper bound on the set of values returned by operator().
  //!
  //!Throws: nothing.
  static result_type max BOOST_PREVENT_MACRO_SUBSTITUTION () { return m - 1u; }

  //!Returns: The dimension of of the quasi-random domain.
  //!
  //!Throws: nothing.
  static std::size_t dimension() { return Dimension; }

  //!Effects: Constructs the default Niederreiter base 2 quasi-random number generator.
  //!
  //!Throws: nothing.
  niederreiter_base2()
    : base_t(std::size_t()) // initialize lattice here
  {}

  //!Effects: Constructs the Niederreiter base 2 quasi-random number generator,
  //!equivalent to X u; u.seed(init).
  //!
  //!Throws: overflow_error.
  explicit niederreiter_base2(std::size_t init)
    : base_t(init) // initialize lattice here
  {}

  //!Requirements: *this is mutable.
  //!
  //!Effects: Resets the quasi-random number generator state to
  //!the one given by the default construction. Equivalent to u.seed(0).
  //!
  //!\brief Throws: nothing.
  void seed()
  {
    base_t::reset_state();
  }

  //!Requirements: *this is mutable.
  //!
  //!Effects: Effectively sets the quasi-random number generator state to the init-th
  //!vector in the s-dimensional quasi-random domain, where s == X::dimension().
  //!\code
  //!X u, v;
  //!for(int i = 0; i < N; ++i)
  //!    for( std::size_t j = 0; j < u.dimension(); ++j )
  //!        u();
  //!v.seed(N);
  //!assert(u() == v());
  //!\endcode
  //!
  //!\brief Throws: overflow_error.
  void seed(std::size_t init)
  {
    base_t::seed(init, "niederreiter_base2::seed");
  }

  //=========================Doxygen needs this!==============================

  //!Requirements: *this is mutable.
  //!
  //!Returns: Returns a successive element of an s-dimensional
  //!(s = X::dimension()) vector at each invocation. When all elements are
  //!exhausted, X::operator() begins anew with the starting element of a
  //!subsequent s-dimensional vector.
  //!
  //!Throws: overflow_error.
  result_type operator()()
  {
    return base_t::operator()();
  }

  //!Requirements: *this is mutable.
  //!
  //!Effects: Advances *this state as if z consecutive
  //!X::operator() invocations were executed.
  //!\code
  //!X u = v;
  //!for(int i = 0; i < N; ++i)
  //!    u();
  //!v.discard(N);
  //!assert(u() == v());
  //!\endcode
  //!
  //!Throws: overflow_error.
  void discard(std::size_t z)
  {
    base_t::discard(z);
  }

  //!Writes a @c niederreiter_base2 to a @c std::ostream.
  BOOST_RANDOM_DETAIL_OSTREAM_OPERATOR(os, niederreiter_base2, f)
  { return f.write_impl(os); }
  //!Reads a @c niederreiter_base2 from a @c std::istream.
  BOOST_RANDOM_DETAIL_ISTREAM_OPERATOR(is, niederreiter_base2, f)
  { return f.read_impl(is); }

  //!Returns true if the two generators will produce identical sequences.
  BOOST_RANDOM_DETAIL_EQUALITY_OPERATOR(niederreiter_base2, x, y)
  { return x.is_equal(y); }
  //!Returns true if the two generators will produce different sequences,
  BOOST_RANDOM_DETAIL_INEQUALITY_OPERATOR(niederreiter_base2)
};

} // namespace random

//It would have been nice to do something like this, but it seems that doxygen
//simply won't show those typedefs, so we spell them out by hand.

/*
#define BOOST_NIEDERREITER_BASE2_TYPEDEF(z, N, text) \
typedef random::niederreiter_base2<uint32_t, N, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_##N##d; \
//

BOOST_PP_REPEAT_FROM_TO(1, 21, BOOST_NIEDERREITER_BASE2_TYPEDEF, unused)

#undef BOOST_NIEDERREITER_BASE2_TYPEDEF
*/

typedef random::niederreiter_base2<uint32_t, 1, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_1d;
typedef random::niederreiter_base2<uint32_t, 2, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_2d;
typedef random::niederreiter_base2<uint32_t, 3, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_3d;
typedef random::niederreiter_base2<uint32_t, 4, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_4d;
typedef random::niederreiter_base2<uint32_t, 5, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_5d;
typedef random::niederreiter_base2<uint32_t, 6, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_6d;
typedef random::niederreiter_base2<uint32_t, 7, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_7d;
typedef random::niederreiter_base2<uint32_t, 8, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_8d;
typedef random::niederreiter_base2<uint32_t, 9, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_9d;
typedef random::niederreiter_base2<uint32_t, 10, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_10d;
typedef random::niederreiter_base2<uint32_t, 11, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_11d;
typedef random::niederreiter_base2<uint32_t, 12, 1, 2*((uint32_t)1 << 31)> niederreiter_base2_12d;

} // namespace boost

#endif // BOOST_RANDOM_NIEDERREITER_BASE2_HPP
