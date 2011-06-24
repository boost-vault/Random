// Copyright Justinas Vygintas Daugmaudis, 2010.
// Use, modification and distribution is subject to the
// Boost Software License, Version 1.0. (See accompanying
// file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)

#include "qrng_typedefs.hpp"

#include <gsl/gsl_qrng.h>


#include <boost/random/uniform_real.hpp>

#include <boost/preprocessor/repeat_from_to.hpp>

#include <boost/shared_ptr.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

//
// DESCRIPTION:
// ~~~~~~~~~~~~
//
// This file tests the quasi-random number generators.
// There are two sets of tests:
// 1) validate_sobol;
// 2) validate_niederreiter_base2.
// These tests compare our results with values produced by qrngs in
// GNU Scientific Library.

// GSL adaptors

struct gsl_niederreiter_base2
{
  static const gsl_qrng_type * type() {
    return ::gsl_qrng_niederreiter_2;
  }
};

struct gsl_sobol
{
  static const gsl_qrng_type * type() {
    return ::gsl_qrng_sobol;
  }
};

struct gsl_deleter
{
  void operator()(gsl_qrng *q) const
  {
    if( q != 0 )
      gsl_qrng_free(q);
  }
};

// Test routine
template<typename GSL, typename QRNG>
bool test_qrng(const std::size_t D)
{
  QRNG test;

  const int validating_sequence_lenght = 100000;

  boost::shared_ptr<gsl_qrng> q(gsl_qrng_alloc(GSL::type(), D), gsl_deleter());
  // RAII, to ensure that no memory leaks occur
  BOOST_REQUIRE( q.get() != 0 );

  double val[D];
  boost::uniform_real<double> rnd;

  for( int i = 0; i < validating_sequence_lenght; ++i )
  {
    gsl_qrng_get(q.get(), val);
    for( std::size_t j = 0; j != D; ++j )
    {
      double n = rnd(test);

      // We want to check that quasi-random number generator values differ no
      // more than 0.0001% of their value.
      BOOST_REQUIRE_CLOSE_FRACTION(val[j], n, 0.0001);
    }
  }

  return true;
}


template<int D>
void test_sobol()
{
  typedef typename sobol_generator<D>::type sobol_t;
  BOOST_CHECK( (test_qrng<gsl_sobol, sobol_t>(D)) );
}

template<int D>
void test_niederreiter_base2()
{
  typedef typename niederreiter_base2_generator<D>::type nb2_t;
  BOOST_CHECK( (test_qrng<gsl_niederreiter_base2, nb2_t>(D)) );
}


BOOST_AUTO_TEST_CASE( validate_sobol )
{
#define UNIT_TEST_QRNG_VALIDATE_SOBOL(z, N, text) test_sobol<N>();

  BOOST_PP_REPEAT_FROM_TO(1, 41, UNIT_TEST_QRNG_VALIDATE_SOBOL, unused)

#undef UNIT_TEST_QRNG_VALIDATE_SOBOL
}

BOOST_AUTO_TEST_CASE( validate_niederreiter_base2 )
{
#define UNIT_TEST_QRNG_VALIDATE_NIEDERREITER_BASE2(z, N, text) test_niederreiter_base2<N>();

  // GSL supports only up to 12 dimensions, so we check only that much
  BOOST_PP_REPEAT_FROM_TO(1, 13, UNIT_TEST_QRNG_VALIDATE_NIEDERREITER_BASE2, unused)

#undef UNIT_TEST_QRNG_VALIDATE_NIEDERREITER_BASE2
}

