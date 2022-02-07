#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/virtualizedgrid.hh>

#include "gridcheck.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  {
    Dune::YaspGrid<1> yaspgrid({1.}, {32});

    Dune::VirtualizedGrid<1, 1> vgrid( yaspgrid );
    // gridcheck(vgrid);

    auto b = vgrid.template leafbegin<0>();
    auto e = vgrid.template leafend<0>();
    std::cout << b->geometry().center() << std::endl;
    // std::cout << e->geometry().center() << std::endl;

    // for( const auto& v : vertices(vgrid.levelGridView(0)) )
    //   std::cout << v.geometry().center() << std::endl;
    //
    // for( const auto& e : elements(vgrid.leafGridView()) )
    //   std::cout << e.geometry().center() << std::endl;
  }

  return 0;
}
