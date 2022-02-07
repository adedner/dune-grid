#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/common/timer.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/virtualizedgrid.hh>

#include "gridcheck.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  {
    Dune::YaspGrid<1> yaspgrid({1.}, {32});
    Dune::VirtualizedGrid<1, 1> vgrid( yaspgrid );

    Dune::Timer timer;
    gridcheck(vgrid);
    std::cout << "Standard grid implementation: " << timer.elapsed() << std::endl;

    timer.reset();
    gridcheck(vgrid);
    std::cout << "Virtualized grid implementation: " << timer.elapsed() << std::endl;
  }

  return 0;
}
