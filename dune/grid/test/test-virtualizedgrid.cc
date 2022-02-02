#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/virtualizedgrid.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  {
    Dune::YaspGrid<2> yaspgrid({1., 1.}, {4, 4});
    Dune::VirtualizedGrid<2, 2> vgrid( yaspgrid );
  }

  // {
  //   Dune::YaspGrid<3> grid({1., 1., 1.}, {4, 4, 4});
  // }

  return 0;
}
