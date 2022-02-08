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


    // // 2D
    // std::cout << "============= 2D =============" << std::endl;
    //
    // Dune::YaspGrid<2> yaspgrid2({1., 1.}, {32, 32});
    // Dune::VirtualizedGrid<2, 2> vgrid2( yaspgrid2 );
    //
    // timer.reset();
    // gridcheck(yaspgrid2);
    // std::cout << "------------------------------" << std::endl;
    // std::cout << "YaspGrid<2>: " << timer.elapsed() << std::endl;
    // std::cout << "------------------------------" << std::endl;
    //
    // timer.reset();
    // gridcheck(vgrid2);
    // std::cout << "------------------------------" << std::endl;
    // std::cout << "Virtualized<2>: " << timer.elapsed() << std::endl;
    // std::cout << "=============================" << std::endl;
    // std::cout << std::endl;
    //
    //
    // // 3D
    // std::cout << "============= 3D =============" << std::endl;
    //
    // Dune::YaspGrid<3> yaspgrid3({1., 1., 1.}, {10, 10, 10});
    // Dune::VirtualizedGrid<3, 3> vgrid3( yaspgrid3 );
    //
    // timer.reset();
    // gridcheck(yaspgrid3);
    // std::cout << "------------------------------" << std::endl;
    // std::cout << "YaspGrid<3>: " << timer.elapsed() << std::endl;
    // std::cout << "------------------------------" << std::endl;
    //
    // timer.reset();
    // gridcheck(vgrid3);
    // std::cout << "------------------------------" << std::endl;
    // std::cout << "Virtualized<3>: " << timer.elapsed() << std::endl;
    // std::cout << "=============================" << std::endl;
  }

  return 0;
}
