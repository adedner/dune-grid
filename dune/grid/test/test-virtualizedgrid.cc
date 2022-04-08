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

  // ======== SMALL TEST OF VIRTUALIZED CLASS ABOVE ========

  int N = 100000;
  Dune::Timer timer;

  // Construct
  timer.reset();
  for (volatile int i = 0; i < N; ++i)
    std::vector<double> v (i, 42.);
  std::cout << "Standard: " << timer.elapsed() << std::endl;

  timer.reset();
  for (volatile int i = 0; i < N; ++i)
    Virtualized virt( std::vector<double> (i, 42.) );
  std::cout << "Virtualized: " << timer.elapsed() << std::endl;

  // Call
  std::vector<double> v (100, 42.);
  timer.reset();
  for (volatile int i = 0; i < N; ++i)
    v.size();
  std::cout << "Call Standard: " << timer.elapsed() << std::endl;

  Virtualized virt(v);
  timer.reset();
  for (volatile int i = 0; i < N; ++i)
    virt.size();
  std::cout << "Call Virtualized: " << timer.elapsed() << std::endl;

  // ===============================================================

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
