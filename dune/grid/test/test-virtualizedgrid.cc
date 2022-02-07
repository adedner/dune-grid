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
    gridcheck(vgrid);


    // auto levelGridView = vgrid.levelGridView(vgrid.maxLevel());
    // auto entity = *levelGridView.template begin<0, Dune::All_Partition>();
    //
    // const auto& globalid = vgrid.globalIdSet();
    //
    // const auto endit = vgrid.leafGridView().iend( entity );
    // auto it = vgrid.leafGridView().ibegin( entity );
    // for(; it != endit; ++it)
    // {
    //   const auto& is = *it;
    //   // std::cout << is.inside().geometry().center() << std::endl;
    //   std::cout << globalid.id(is.inside()) << std::endl;
    // }

    //
    // auto seed = entity.seed();
    //
    //
    // using ThisType = Dune::VirtualizedGrid<1, 1>;
    // using I = Dune::YaspGrid<1>;
    // typedef typename Dune::VirtualizedGridEntitySeed<0, const ThisType>::template Implementation<typename std::decay_t<I>::template Codim<0>::EntitySeed> ImplSeed0;
    // // auto e = yaspgrid.entity(
    // //   dynamic_cast<const ImplSeed0*>(seed.impl().impl_.get())->impl()
    // // );
    // // auto e = vgrid.entity( seed );
    //
    //
    // try {
    //   auto e = dynamic_cast<ImplSeed0&>(*seed.impl().impl_.get());
    // } catch(const std::bad_cast& e)
    // {
    //     std::cout << e.what() << '\n';
    // }

    // seed.impl().impl_.get())->impl()
    //
    // // vgrid.mark( 1, entity );
    //
    // auto it = vgrid.levelGridView(vgrid.maxLevel()).template begin<0>();
    // const auto endit = vgrid.levelGridView(vgrid.maxLevel()).template end<0>();

    // for (; it!=endit; ++it)
    // {
    //   std::cout << it->geometry().center() << std::endl;
    //   std::cout << levelIndexSet.index( *it ) << std::endl;
    //
    //   auto subEntity = it->template subEntity< 0 >( 0 );
    //   std::cout << subEntity.geometry().center() << std::endl;
    //   std::cout << levelIndexSet.index( subEntity ) << std::endl;
    // }

    // using ImplEntity = Dune::YaspGrid<1>::template Codim<0>::Entity;
    // using I = Dune::YaspGrid<1>::template Codim<0>::Entity;
    // using V = Dune::VirtualizedGrid<1, 1>;
    // typedef typename Dune::VirtualizedGridEntity<0, V::dimension, V>::Implementation<I> ImplEntity;
    // const ImplEntity* eImpl = dynamic_cast<const ImplEntity*>(entity.impl().impl_.get())->impl();

    // std::cout << vgrid.leafIndexSet().index(entity) << std::endl;

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
