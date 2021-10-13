#include <config.h>

#include <iostream>
#include <set>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/entityowner.hh>

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
}


// Check whether the extracted ownership is correct
void testOwnership (Dune::TestSuite& testSuite,
                    int rank, int ownerRank,
                    Dune::PartitionType pt, Dune::OwnerType ot)
{
  if (pt == Dune::PartitionType::InteriorEntity) {
    testSuite.check(ownerRank == rank);
    testSuite.check(ot == Dune::OwnerType::Owner);
  } else if (pt == Dune::PartitionType::OverlapEntity) {
    testSuite.check(ownerRank != rank);
    testSuite.check(ot == Dune::OwnerType::Overlap);
  } else if(pt == Dune::PartitionType::FrontEntity || pt == Dune::PartitionType::GhostEntity)  {
    testSuite.check(ownerRank != rank);
    testSuite.check(ot == Dune::OwnerType::Ghost);
  } else if (pt == Dune::PartitionType::BorderEntity) {
    testSuite.check(ot == Dune::OwnerType::Owner || ot == Dune::OwnerType::Overlap);
    // TODO: check that exactly one process is owner and all others are overlap
  }
}


// check the elements of a given gridView for correct ownership computation
template <class GridView>
void test (std::string name, Dune::TestSuite& testSuite, GridView const& gridView)
{
  Dune::TestSuite subTestSuite{name};

  print(name, ":\n");
  const int rank = gridView.comm().rank();

  Dune::EntityOwner entityOwner{gridView};
  for (auto const& e : elements(gridView))
  {
    testOwnership(subTestSuite, rank, entityOwner.owner(e),
                  e.partitionType(), entityOwner.ownerType(e));

    Dune::Hybrid::forEach(Dune::StaticIntegralRange<int,GridView::dimension+1>{}, [&](auto codim)
    {
      for (unsigned int i = 0; i < e.subEntities(codim); ++i)
      {
        Dune::PartitionType pt = e.template subEntity<codim>(i).partitionType();
        int ownerRank = entityOwner.owner(e,i,codim);
        Dune::OwnerType ot = entityOwner.ownerType(e,i,codim);

        testOwnership(subTestSuite, rank, ownerRank, pt, ot);
      }
    });
  }

  testSuite.subTest(subTestSuite);
}


int main (int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite{"EntityOwnerTest"};

  { // Test YaspGrid 2d
    Dune::YaspGrid<2> grid2d0o({1.0,1.0}, {4,4}, 0u, 0); // overlap = 0
    test("YaspGrid<2> overlap=0", testSuite, grid2d0o.leafGridView());
    Dune::YaspGrid<2> grid2d1o({1.0,1.0}, {8,8}, 0u, 1); // overlap = 1
    test("YaspGrid<2> overlap=1", testSuite, grid2d1o.leafGridView());
    Dune::YaspGrid<2> grid2d2o({1.0,1.0}, {16,16}, 0u, 2); // overlap = 2
    test("YaspGrid<2> overlap=2", testSuite, grid2d2o.leafGridView());
  }

  { // Test YaspGrid 3d
    Dune::YaspGrid<3> grid3d0o({1.0,1.0,1.0}, {8,8,8}, 0u, 0); // overlap = 0
    test("YaspGrid<3> overlap=0", testSuite, grid3d0o.leafGridView());
    Dune::YaspGrid<3> grid3d1o({1.0,1.0,1.0}, {16,16,16}, 0u, 1); // overlap = 1
    test("YaspGrid<3> overlap=1", testSuite, grid3d1o.leafGridView());
    Dune::YaspGrid<3> grid3d2o({1.0,1.0,1.0}, {32,32,32}, 0u, 2); // overlap = 2
    test("YaspGrid<3> overlap=2", testSuite, grid3d2o.leafGridView());
  }

#if HAVE_DUNE_UGGRID
  { // Test UGGrid 2d and 3d
    using Factory2d = Dune::StructuredGridFactory<Dune::UGGrid<2>>;
    using Factory3d = Dune::StructuredGridFactory<Dune::UGGrid<3>>;

    auto uggrid2ds = Factory2d::createSimplexGrid({0.0,0.0},{1.0,1.0},{4u,4u});
    uggrid2ds->loadBalance();
    test("UGGrid<2> simplex", testSuite, uggrid2ds->leafGridView());

    auto uggrid2dc = Factory2d::createCubeGrid({0.0,0.0},{1.0,1.0},{4u,4u});
    uggrid2dc->loadBalance();
    test("UGGrid<2> cube", testSuite, uggrid2dc->leafGridView());

    auto uggrid3ds = Factory3d::createSimplexGrid({0.0,0.0,0.0},{1.0,1.0,1.0},{4u,4u,4u});
    uggrid3ds->loadBalance();
    test("UGGrid<3> simplex", testSuite, uggrid3ds->leafGridView());

    auto uggrid3dc = Factory3d::createCubeGrid({0.0,0.0,0.0},{1.0,1.0,1.0},{4u,4u,4u});
    uggrid3dc->loadBalance();
    test("UGGrid<3> cube", testSuite, uggrid3dc->leafGridView());
  }
#endif

  return testSuite.exit();
}