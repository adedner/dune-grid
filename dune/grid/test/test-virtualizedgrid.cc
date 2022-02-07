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

    auto levelGridView = vgrid.levelGridView(0);
    auto entity = *levelGridView.template begin<0, Dune::All_Partition>();
    std::cout << entity.geometry().center() << std::endl;

    // using ImplEntity = Dune::YaspGrid<1>::template Codim<0>::Entity;
    // using I = Dune::YaspGrid<1>::template Codim<0>::Entity;
    // using V = Dune::VirtualizedGrid<1, 1>;
    // typedef typename Dune::VirtualizedGridEntity<0, V::dimension, V>::Implementation<I> ImplEntity;
    // const ImplEntity* eImpl = dynamic_cast<const ImplEntity*>(entity.impl().impl_.get())->impl();

    std::cout << vgrid.leafIndexSet().index(entity) << std::endl;

    // auto i = *levelGridView.template begin<0,Dune::InteriorBorder_Partition>();

    // const auto& idxSet = levelGridView.indexSet();

    // auto idxa = levelGridView.indexSet().index(a);
    // auto idxi = levelGridView.indexSet().index(i);


    // auto b = vgrid.template leafbegin<0>();
    // auto e = vgrid.template leafend<0>();
    // for ( const auto& e : elements(vgrid.leafGridView()) )
      // std::cout << e.geometry().center() << std::endl;
    // std::cout << e->geometry().center() << std::endl;

    // for( const auto& v : vertices(vgrid.levelGridView(0)) )
    //   std::cout << v.geometry().center() << std::endl;
    //
    // for( const auto& e : elements(vgrid.leafGridView()) )
    //   std::cout << e.geometry().center() << std::endl;
  }

  return 0;
}
