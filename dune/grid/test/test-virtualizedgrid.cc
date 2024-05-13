// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/common/timer.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/virtualizedgrid.hh>

#include "gridcheck.hh"

class Virtualized
{
  struct Interface
  {
    virtual ~Interface () = default;
    virtual int size () const = 0;
  };

  template< class I >
  struct Implementation final
    : public Interface
  {
    Implementation (I&& i) : impl_( std::forward<I>(i) ) {}
    virtual int size () const override { return impl_.size(); }

  private:
    I impl_;
  };

public:
  template< class Impl >
  Virtualized(Impl&& impl) : impl_( new Implementation<Impl>( std::forward<Impl>(impl) ) ) {}
  int size () const { return impl_->size(); }
  std::unique_ptr<Interface> impl_;
};

template<class Grid>
void run(std::string name, Grid& grid)
{
  const int dim = Grid::dimension;
  const int dow = Grid::dimensionworld;
  std::cout << "============= " << dim << "D/" << dow << "D =============" << std::endl;

  Dune::VirtualizedGrid<dim,dow> vgrid( grid );

  Dune::Timer timer;
  gridcheck(grid);
  std::cout << "------------------------------" << std::endl;
  std::cout << name << ": " << timer.elapsed()  << std::endl;
  std::cout << "------------------------------" << std::endl;

  timer.reset();
  gridcheck(vgrid);
  std::cout << "------------------------------" << std::endl;
  std::cout << "Virtualized<" << dim << "," << dow << ">: " << timer.elapsed() << std::endl;
  std::cout << "=============================" << std::endl;
  std::cout << std::endl;
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance( argc, argv );

  // ======== SMALL TEST OF VIRTUALIZED CLASS ABOVE ========

  volatile int N = 1000;
  Dune::Timer timer;

  using T = std::vector<double>;

  // Construct
  timer.reset();
  std::vector<T> standard;
  for (int i = 0; i < N; ++i)
    standard.emplace_back( i, 42. );
  std::cout << "Standard: " << timer.elapsed() << std::endl;

  timer.reset();
  std::vector<Virtualized> virtualized;
  for (int i = 0; i < N; ++i)
    virtualized.emplace_back( T(i, 42.) );
  std::cout << "Virtual:  " << timer.elapsed() << std::endl;

  // Call
  timer.reset();
  int c = 0;
  for (int i = 0; i < 10000; ++i)
    for (int i = 0; i < N; ++i)
      c += standard[i].size();
  std::cout << "Call Standard (c = " << c << "): " << timer.elapsed() << std::endl;

  timer.reset();
  c = 0;
  for (int i = 0; i < 10000; ++i)
    for (int i = 0; i < N; ++i)
      c += virtualized[i].size();
  std::cout << "Call Virtual  (c = " << c << "): " << timer.elapsed() << std::endl;

  // ======== SMALL TEST COMPUTING CUMULATED VOLUME ========

  Dune::YaspGrid<2> yaspgrid2({1., 1.}, {1000, 1000});
  Dune::VirtualizedGrid<2, 2> vgrid2( yaspgrid2 );

  // Compute volume
  std::cout << "Volume calculation" << std::endl;
  timer.reset();
  volatile double vol = 0.0;
  for (const auto& e : elements( yaspgrid2.leafGridView() ))
    vol += e.geometry().volume();
  std::cout << " Standard (vol = " << vol << "): " << timer.elapsed() << std::endl;

  timer.reset();
  vol = 0.0;
  for (const auto& e : elements( vgrid2.leafGridView() ))
    vol += e.geometry().volume();
  std::cout << " Virtual  (vol = " << vol << "): " << timer.elapsed() << std::endl;

  // ======== GRID CHECK ========

  {
    Dune::OneDGrid grid1(32, 0.0, 1.0);
    run("OneDGrid", grid1);
  }

  {
    Dune::YaspGrid<1> yaspgrid1({1.}, {32});
    run("YaspGrid<1>", yaspgrid1);

    Dune::YaspGrid<2> yaspgrid2({1., 1.}, {6, 6});
    run("YaspGrid<2>", yaspgrid2);

    Dune::YaspGrid<3> yaspgrid3({1., 1., 1.}, {4, 4, 4});
    run("YaspGrid<3>", yaspgrid3);
  }

#if HAVE_ALBERTA
  {
    using Grid1 = Dune::AlbertaGrid<1,3>;
    using Factory1 = Dune::StructuredGridFactory<Grid1>;
    auto grid1ptr = Factory1::createSimplexGrid({0.0,0.0,0.0}, {1.0,0.0,0.0}, {32u});
    run("AlbertaGrid<1,3>", *grid1ptr);

    using Grid2 = Dune::AlbertaGrid<2,3>;
    using Factory2 = Dune::StructuredGridFactory<Grid2>;
    auto grid2ptr = Factory2::createSimplexGrid({0.0,0.0,0.0}, {1.0,1.0,0.0}, {32u,32u});
    run("AlbertaGrid<2,3>", *grid2ptr);

    using Grid3 = Dune::AlbertaGrid<3,3>;
    using Factory3 = Dune::StructuredGridFactory<Grid3>;
    auto grid3ptr = Factory3::createSimplexGrid({0.0,0.0,0.0}, {1.0,1.0,1.0}, {32u,32u,32u});
    run("AlbertaGrid<3,3>", *grid3ptr);
  }
#endif // HAVE_ALBERTA

#if HAVE_DUNE_UGGRID
  {
    using Grid2 = Dune::UGGrid<2>;
    using Factory2 = Dune::StructuredGridFactory<Grid2>;
    auto grid2ptr = Factory2::createSimplexGrid({0.0,0.0}, {1.0,1.0}, {32u,32u});
    run("UGGrid<2>", *grid2ptr);

    using Grid3 = Dune::UGGrid<3>;
    using Factory3 = Dune::StructuredGridFactory<Grid3>;
    auto grid3ptr = Factory3::createSimplexGrid({0.0,0.0,0.0}, {1.0,1.0,1.0}, {32u,32u,32u});
    run("UGGrid<3>", *grid3ptr);
  }
#endif

  return 0;
}
