#include <config.h>

#include <iostream>
#include <set>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/entityowner.hh>
#include <dune/grid/utility/entityownertype.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

// some helper functions for parallel output
namespace
{
  template <class Arg, class... Args>
  std::string concat (Arg&& arg, Args&&... args)
  {
    std::stringstream ss;
    ss << std::forward<Arg>(arg);
    if constexpr(sizeof...(args) > 0)
      ss << concat(std::forward<Args>(args)...);
    return ss.str();
  }

  template <class... Args>
  void print (Args&&... args)
  {
    Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance();
    std::cout << concat("[", mpiHelper.rank(), "]  ", std::forward<Args>(args)...) << std::flush;
  }

  template <class... Args>
  void println (Args&&... args)
  {
    Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance();
    std::cout << concat("[", mpiHelper.rank(), "]  ", std::forward<Args>(args)...) << std::endl;
  }
}

class EntityOwnerTest
{
public:
  EntityOwnerTest (const Dune::MPIHelper& mpiHelper)
    : testSuite_{"EntityOwner"}
    , rank_{mpiHelper.rank()}
    , size_{mpiHelper.size()}
  {}

  template <class GridView>
  void test (std::string name, const GridView& gridView)
  {
    println(name, ":");

    Dune::EntityOwner entityOwner{gridView};
    Dune::EntityOwnerType entityOwnerType{gridView};
    for (auto const& e : elements(gridView))
    {
      testOwnership(entityOwner.ownerRank(e), e.partitionType(), entityOwner.ownerType(e));
      testSuite_.check(entityOwner.ownerType(e) == entityOwnerType.ownerType(e));

      Dune::Hybrid::forEach(Dune::StaticIntegralRange<int,GridView::dimension+1>{}, [&](auto codim)
      {
        for (unsigned int i = 0; i < e.subEntities(codim); ++i)
        {
          Dune::PartitionType pt = e.template subEntity<codim>(i).partitionType();
          int ownerRank = entityOwner.ownerRank(e,i,codim);
          Dune::OwnerType ot = entityOwner.ownerType(e,i,codim);

          testOwnership(ownerRank, pt, ot);
          testSuite_.check(ot == entityOwnerType.ownerType(e,i,codim));
        }
      });
    }
  }

  // Check whether the extracted ownership is correct
  void testOwnership (int ownerRank, Dune::PartitionType pt, Dune::OwnerType ot)
  {
    if (pt == Dune::PartitionType::InteriorEntity) {
      testSuite_.check(ownerRank == rank_, "InteriorEntity");
      testSuite_.check(ot == Dune::OwnerType::Owner, "InteriorEntity");
    } else if (pt == Dune::PartitionType::OverlapEntity) {
      testSuite_.check(ownerRank != rank_, "OverlapEntity");
      testSuite_.check(ot == Dune::OwnerType::Overlap, "OverlapEntity");
    } else if(pt == Dune::PartitionType::FrontEntity || pt == Dune::PartitionType::GhostEntity)  {
      testSuite_.check(ownerRank != rank_, "GhostEntity");
      testSuite_.check(ot == Dune::OwnerType::Ghost, "GhostEntity");
    } else if (pt == Dune::PartitionType::BorderEntity) {
      testSuite_.check(ot == Dune::OwnerType::Owner ||
                       ot == Dune::OwnerType::Overlap, "BorderEntity");
    }
  }

  int exit () const
  {
    return testSuite_.exit();
  }

private:
  Dune::TestSuite testSuite_;
  int const rank_;
  int const size_;
};


int main (int argc, char** argv)
{
  const Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);

  EntityOwnerTest testSuite{mpiHelper};

  { // Test YaspGrid 2d
    Dune::YaspGrid<2> grid2d0o({1.0,1.0}, {4,4}, 0u, 0); // overlap = 0
    testSuite.test("YaspGrid<2> overlap=0", grid2d0o.leafGridView());
    Dune::YaspGrid<2> grid2d1o({1.0,1.0}, {8,8}, 0u, 1); // overlap = 1
    testSuite.test("YaspGrid<2> overlap=1", grid2d1o.leafGridView());
    Dune::YaspGrid<2> grid2d2o({1.0,1.0}, {16,16}, 0u, 2); // overlap = 2
    testSuite.test("YaspGrid<2> overlap=2", grid2d2o.leafGridView());
  }

  { // Test YaspGrid 3d
    Dune::YaspGrid<3> grid3d0o({1.0,1.0,1.0}, {8,8,8}, 0u, 0); // overlap = 0
    testSuite.test("YaspGrid<3> overlap=0", grid3d0o.leafGridView());
    Dune::YaspGrid<3> grid3d1o({1.0,1.0,1.0}, {16,16,16}, 0u, 1); // overlap = 1
    testSuite.test("YaspGrid<3> overlap=1", grid3d1o.leafGridView());
    Dune::YaspGrid<3> grid3d2o({1.0,1.0,1.0}, {32,32,32}, 0u, 2); // overlap = 2
    testSuite.test("YaspGrid<3> overlap=2", grid3d2o.leafGridView());
  }

#if HAVE_DUNE_UGGRID
  { // Test UGGrid 2d and 3d
    using Factory2d = Dune::StructuredGridFactory<Dune::UGGrid<2>>;
    using Factory3d = Dune::StructuredGridFactory<Dune::UGGrid<3>>;

    auto uggrid2ds = Factory2d::createSimplexGrid({0.0,0.0},{1.0,1.0},{4u,4u});
    uggrid2ds->loadBalance();
    testSuite.test("UGGrid<2> simplex", uggrid2ds->leafGridView());

    auto uggrid2dc = Factory2d::createCubeGrid({0.0,0.0},{1.0,1.0},{4u,4u});
    uggrid2dc->loadBalance();
    testSuite.test("UGGrid<2> cube", uggrid2dc->leafGridView());

    auto uggrid3ds = Factory3d::createSimplexGrid({0.0,0.0,0.0},{1.0,1.0,1.0},{4u,4u,4u});
    uggrid3ds->loadBalance();
    testSuite.test("UGGrid<3> simplex", uggrid3ds->leafGridView());

    auto uggrid3dc = Factory3d::createCubeGrid({0.0,0.0,0.0},{1.0,1.0,1.0},{4u,4u,4u});
    uggrid3dc->loadBalance();
    testSuite.test("UGGrid<3> cube", uggrid3dc->leafGridView());
  }
#endif

  return testSuite.exit();
}