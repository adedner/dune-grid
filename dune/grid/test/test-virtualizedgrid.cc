// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/common/timer.hh>
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
    // 1D
    std::cout << "============= 1D =============" << std::endl;

    Dune::YaspGrid<1> yaspgrid1({1.}, {32});
    Dune::VirtualizedGrid<1, 1> vgrid1( yaspgrid1 );

    Dune::Timer timer;
    gridcheck(yaspgrid1);
    std::cout << "------------------------------" << std::endl;
    std::cout << "YaspGrid<1>: " << timer.elapsed() << std::endl;
    std::cout << "------------------------------" << std::endl;

    timer.reset();
    gridcheck(vgrid1);
    std::cout << "------------------------------" << std::endl;
    std::cout << "Virtualized<1>: " << timer.elapsed() << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;


    // 2D
    std::cout << "============= 2D =============" << std::endl;

    Dune::YaspGrid<2> yaspgrid2({1., 1.}, {6, 6});
    Dune::VirtualizedGrid<2, 2> vgrid2( yaspgrid2 );

    timer.reset();
    gridcheck(yaspgrid2);
    std::cout << "------------------------------" << std::endl;
    std::cout << "YaspGrid<2>: " << timer.elapsed() << std::endl;
    std::cout << "------------------------------" << std::endl;

    timer.reset();
    gridcheck(vgrid2);
    std::cout << "------------------------------" << std::endl;
    std::cout << "Virtualized<2>: " << timer.elapsed() << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;


    // 3D
    std::cout << "============= 3D =============" << std::endl;

    Dune::YaspGrid<3> yaspgrid3({1., 1., 1.}, {4, 4, 4});
    Dune::VirtualizedGrid<3, 3> vgrid3( yaspgrid3 );

    timer.reset();
    gridcheck(yaspgrid3);
    std::cout << "------------------------------" << std::endl;
    std::cout << "YaspGrid<3>: " << timer.elapsed() << std::endl;
    std::cout << "------------------------------" << std::endl;

    timer.reset();
    gridcheck(vgrid3);
    std::cout << "------------------------------" << std::endl;
    std::cout << "Virtualized<3>: " << timer.elapsed() << std::endl;
    std::cout << "=============================" << std::endl;
  }

  return 0;
}
