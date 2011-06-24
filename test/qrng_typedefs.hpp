// Copyright Justinas Vygintas Daugmaudis, 2010.
// Use, modification and distribution is subject to the
// Boost Software License, Version 1.0. (See accompanying
// file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)

#ifndef HEADER_QRNG_TYPEDEFS_HPP_INCLUDED
#define HEADER_QRNG_TYPEDEFS_HPP_INCLUDED

#include <boost/random/niederreiter_base2.hpp>
#include <boost/random/sobol.hpp>
#include <boost/random/faure.hpp>


// Boost.Random adaptors

template<int D>
struct niederreiter_base2_generator
{
  typedef boost::random::niederreiter_base2<boost::uint32_t, D, 1, 2*((boost::uint32_t)1 << 31)> type;
};

template<int D>
struct sobol_generator
{
  typedef boost::random::sobol<boost::uint32_t, D, 1, (boost::uint32_t)1 << 31> type;
};

template<int D>
struct faure_generator
{
  typedef boost::random::faure<double, D> type;
};


// To be used in place of uniform_real, where scaling is not necessary.
// Identity distribution simply returns the values as produced by an
// engine.
template<typename T>
struct identity_distribution
{
  typedef T result_type;
  template<typename QRNG>
  result_type operator()(QRNG& q) const
  { return q(); }
};

#endif // HEADER_QRNG_TYPEDEFS_HPP_INCLUDED
