// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/* Tests for periodic YaspGrids (and the overlap sizes they are used with)
 *
 * A periodic YaspGrid stores copies of the entities at the opposite side of the
 * domain in its overlap, also on a single process. The copies have unwrapped
 * geometries (outside of the domain), their own indices, but the same ids as
 * the originals. The grid checks are told about this by a periodic map, see
 * dune/grid/test/periodicmap.hh.
 */

#include <bitset>
#include <iostream>
#include <string>

#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include "test-yaspgrid.hh"

// all coordinate types, periodic in direction 0, overlap 1, with and without refinement
template<int dim>
void checkCoordinateTypes (const std::string& testID)
{
  using Equidistant = Dune::EquidistantCoordinates<double,dim>;
  using EquidistantOffset = Dune::EquidistantOffsetCoordinates<double,dim>;
  using TensorProduct = Dune::TensorProductCoordinates<double,dim>;

  for (int refine = 0; refine <= 1; ++refine)
  {
    const std::string id = testID + "-ref" + std::to_string(refine);
    check_yasp(id + "equidistant",
               YaspFactory<dim,Equidistant>::buildGrid(true, refine, true),
               yaspPeriodicMap<dim,Equidistant>(true));
    check_yasp(id + "equidistantoffset",
               YaspFactory<dim,EquidistantOffset>::buildGrid(true, refine, true),
               yaspPeriodicMap<dim,EquidistantOffset>(true));
    check_yasp(id + "tensor",
               YaspFactory<dim,TensorProduct>::buildGrid(true, refine, true),
               yaspPeriodicMap<dim,TensorProduct>(true));
  }
}

// overlap 0, 1, 2 for no, one, and all periodic directions
template<int dim>
void checkOverlapAndDirections (const std::string& testID, int ranks)
{
  using Coordinates = Dune::EquidistantOffsetCoordinates<double,dim>;
  using Grid = Dune::YaspGrid<dim,Coordinates>;

  const Dune::FieldVector<double,dim> lower(-1.0), upper(1.0);
  std::array<int,dim> size;
  std::fill(size.begin(), size.end(), 4);
  size[0] = 8;

  const std::bitset<dim> none, first(1), all((1 << dim) - 1);

  // overlap 2 does not fit into this grid with more than 2 processes
  const int maxOverlap = (ranks > 2) ? 1 : 2;

  for (int overlap = 0; overlap <= maxOverlap; ++overlap)
    for (const auto& periodic : { none, first, all })
      for (int refine = 0; refine <= 1; ++refine)
      {
        const std::string id = testID + "-overlap" + std::to_string(overlap)
          + "-periodic" + periodic.to_string() + "-ref" + std::to_string(refine);
        std::cout << std::endl << "Checking " << id << std::endl;

        Grid* grid = new Grid(lower, upper, size, periodic, overlap);
        grid->globalRefine(refine);
        check_yasp(id, grid, Dune::GridCheck::PeriodicMap<Grid>(
                     Dune::GridCheck::periodicBoxMap(lower, upper, periodic)));
      }
}

// coarse grids: with two cells in a periodic direction, different cells have the
// same set of mapped corners, so the grid checks must not identify entities by
// their corners alone
template<int dim>
void checkCoarse (const std::string& testID, int ranks)
{
  if (ranks > 1)
    return;

  using Coordinates = Dune::EquidistantOffsetCoordinates<double,dim>;
  using Grid = Dune::YaspGrid<dim,Coordinates>;

  const Dune::FieldVector<double,dim> lower(-1.0), upper(1.0);
  std::array<int,dim> size;
  std::fill(size.begin(), size.end(), 2);

  const std::bitset<dim> first(1), all((1 << dim) - 1);
  for (const auto& periodic : { first, all })
  {
    const std::string id = testID + "-coarse-periodic" + periodic.to_string();
    std::cout << std::endl << "Checking " << id << std::endl;

    Grid* grid = new Grid(lower, upper, size, periodic, 1);
    check_yasp(id, grid, Dune::GridCheck::PeriodicMap<Grid>(
                 Dune::GridCheck::periodicBoxMap(lower, upper, periodic)));
  }
}

int main (int argc , char **argv) {
  try {
    // Initialize MPI, if present
    const auto & mpiHelper = Dune::MPIHelper::instance(argc, argv);

    Dune::Hybrid::forEach(std::integer_sequence<int,1,2,3>{}, [&](auto dim) {
      const std::string testID = "yasp-periodic-" + std::to_string(dim) + "d-np" + std::to_string(mpiHelper.size());
      checkCoordinateTypes<dim>(testID);
      checkOverlapAndDirections<dim>(testID, mpiHelper.size());
      checkCoarse<dim>(testID, mpiHelper.size());
    });

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
