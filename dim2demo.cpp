// Copyright Justinas Vygintas Daugmaudis, 2010.
// Use, modification and distribution is subject to the
// Boost Software License, Version 1.0. (See accompanying
// file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)


#include <boost/random/sobol.hpp>
#include <boost/random/niederreiter_base2.hpp>
#include <boost/random/faure.hpp>

#include <boost/random/uniform_real.hpp>

#include <boost/gil/image.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/gil/extension/io/png_io.hpp>

using namespace boost::gil;

// Models a Unary Function
template <typename Qrng, typename P>   // Models PixelValueConcept
struct pixelnoise_fn {
  typedef point2<ptrdiff_t>    point_t;

  typedef pixelnoise_fn        const_t;
  typedef P                    value_type;
  typedef value_type           reference;
  typedef value_type           const_reference;
  typedef point_t              argument_type;
  typedef reference            result_type;
  BOOST_STATIC_CONSTANT(bool, is_mutable=false);

  Qrng gen;
  boost::uniform_real<> rnd;

  value_type                    _pnt_color;
  point_t                       _img_size;
  std::vector<point_t>          _random_points;

  static const int MAX_ITER=6000;        // max number of iterations

  pixelnoise_fn() {}
  pixelnoise_fn(const point_t& sz, const value_type& in_color) : _pnt_color(in_color), _img_size(sz) { gen_random(); }

  result_type operator()(const point_t& p) const
  {
    // Admitedly silly and slow, but this is straightforward and enough for demonstration purposes
    return std::find( _random_points.begin(), _random_points.end(), p ) != _random_points.end() ? _pnt_color: value_type(255,255,255);
  }

private:
  void gen_random()
  {
    for( int i = 0; i < MAX_ITER; ++i )
      _random_points.push_back( point_t( _img_size.x * rnd(gen), _img_size.y * rnd(gen) ) );
  }
};


template<typename Qrng>
void output_image(const char* img_name)
{
  typedef pixelnoise_fn<Qrng, rgb8_pixel_t> deref_t;
  typedef typename deref_t::point_t  point_t;
  typedef virtual_2d_locator<deref_t,false> locator_t;
  typedef image_view<locator_t> my_virt_view_t;

  boost::function_requires<PixelLocatorConcept<locator_t> >();
  gil_function_requires<StepIteratorConcept<typename locator_t::x_iterator> >();

  point_t dims(500,500);
  my_virt_view_t virt_img(dims, locator_t(point_t(0,0), point_t(1,1), deref_t(dims, rgb8_pixel_t(0,0,0))));
  png_write_view(img_name, virt_img);
}

int main()
{
  output_image<boost::sobol_4d>("out-sobol.png");
  output_image<boost::niederreiter_base2_4d>("out-niederreiter_base2.png");
  output_image<boost::faure_4d>("out-faure.png");
  return 0;
}

