/* boost random/detail/gray_coded_qrng_base.hpp header file
 *
 * Copyright Justinas Vygintas Daugmaudis 2010
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BOOST_RANDOM_DETAIL_GRAY_CODED_QRNG_BASE_HPP
#define BOOST_RANDOM_DETAIL_GRAY_CODED_QRNG_BASE_HPP

#include <boost/random/detail/quasi_random_number_generator_base.hpp>

//!\file
//!Describes the gray-coded quasi-random number generator base class template.

namespace boost {
namespace random {

namespace detail {

template<typename DerivedT, typename LatticeT>
class gray_coded_qrng_base : public quasi_random_number_generator_base<
                                       DerivedT
                                     , const LatticeT
                                    >
{
private:

  typedef gray_coded_qrng_base<DerivedT, LatticeT> self_t;
  typedef quasi_random_number_generator_base<DerivedT, const LatticeT> base_t;

  friend class quasi_random_number_generator_base<
                 DerivedT
               , const LatticeT>; // immutable lattice!

public:

  typedef typename LatticeT::result_type result_type;

  gray_coded_qrng_base()
    : base_t(std::size_t())
  {}

  explicit gray_coded_qrng_base(std::size_t init)
    : base_t(init)
  {}

  // default copy c-tor is fine

  //!Requirements: *this is mutable.
  //!
  //!Effects: *this and rhs are equivalent as they both return
  //!an identical sequence of quasi-random numbers.
  //!\code
  //!X u; u = rhs;
  //!assert(u() == rhs());
  //!\endcode
  //!
  //!Returns: *this.
  //!
  //!Throws: nothing.
  self_t& operator=(const self_t& rhs)
  {
    this->curr_elem = rhs.curr_elem;
    this->seq_count = rhs.seq_count;
    std::copy(rhs.quasi_state, rhs.quasi_state + base_t::dimension_value, this->quasi_state);
    // we do not copy lattice here!
    return *this;
  }

protected:
  void reset_state()
  {
    this->curr_elem = 0;
    this->seq_count = 0;
    std::fill(this->quasi_state, this->quasi_state + base_t::dimension_value, 0);
  }

  void seed(std::size_t init, const char *msg)
  {
    this->curr_elem = 0;
    if ( init != this->seq_count )
    {
      base_t::derived().seed();

      this->seq_count = init;
      init ^= (init / 2);
      for( int r = 0; init != 0; ++r, init >>= 1 )
      {
        if( init & 1 )
          update_quasi(r, msg);
      }
    }
  }

  void compute_current_vector(const char *msg)
  {
    update_state(this->seq_count, msg);
  }

private:

  void compute_next_vector()
  {
    update_state(LatticeT::checked_increment(this->seq_count), "compute_next_vector");
  }

  void update_state(std::size_t cnt, const char *msg)
  {
    // Find the position of the least-significant zero in sequence count.
    // This is the bit that changes in the Gray-code representation as
    // the count is advanced.
    int r = 0;
    for( ; cnt & 1; cnt >>= 1 ) {
      ++r;
    }
    update_quasi(r, msg);
  }

  void update_quasi(int r, const char* msg)
  {
    if(r >= LatticeT::bit_count)
      boost::throw_exception( std::overflow_error(msg) );

    // Calculate the next state.
    for(std::size_t i = 0; i < base_t::dimension_value; ++i)
      this->quasi_state[i] ^= this->lattice(r, i);
  }
};

}} // namespace detail::random

} // namespace boost

#endif // BOOST_RANDOM_DETAIL_GRAY_CODED_QRNG_BASE_HPP
