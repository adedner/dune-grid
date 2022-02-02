#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/virtualizedgrid.hh>

#include "gridcheck.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  {
    Dune::YaspGrid<2> yaspgrid({1., 1.}, {6, 6});

    Dune::VirtualizedGrid<2, 2> vgrid( yaspgrid );
    gridcheck(vgrid);

    // for( const auto& e : elements(vgrid.leafGridView()) )
    //   std::cout << e.geometry().center() << std::endl;
  }

  return 0;
}
