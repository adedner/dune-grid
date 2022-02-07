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
    // 1D
    std::cout << "============= 1D =============" << std::endl;

    Dune::YaspGrid<1> yaspgrid({1.}, {32});
    Dune::VirtualizedGrid<1, 1> vgrid( yaspgrid );

    Dune::Timer timer;
    gridcheck(vgrid);
    std::cout << "------------------------------" << std::endl;
    std::cout << "YaspGrid<1>: " << timer.elapsed() << std::endl;
    std::cout << "------------------------------" << std::endl;

    timer.reset();
    gridcheck(vgrid);
    std::cout << "------------------------------" << std::endl;
    std::cout << "Virtualized<1>: " << timer.elapsed() << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;

    // 2D
    std::cout << "============= 2D =============" << std::endl;

    Dune::YaspGrid<2> yaspgrid2({1., 1.}, {8, 8});
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
  }

  return 0;
}
