// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKIDSET_HH
#define DUNE_GRID_TEST_CHECKIDSET_HH

#include <iostream>
#include <map>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/exceptions.hh>

#include <dune/grid/test/coordinatetolerance.hh>
#include <dune/grid/test/periodicmap.hh>

/** @file
   @brief Unit tests for IdSet implementations
 */


namespace Dune
{
  namespace GridCheck::Impl
  {
    // order (codim, geometry type, center) tuples, comparing centers up to an absolute tolerance
    struct PositionLess
    {
      double tolerance;

      template<class Tuple>
      bool operator() (const Tuple& a, const Tuple& b) const
      {
        if (std::get<0>(a) != std::get<0>(b))
          return std::get<0>(a) < std::get<0>(b);
        if (std::get<1>(a) != std::get<1>(b))
          return std::get<1>(a) < std::get<1>(b);
        return ApproxLess{ tolerance }(std::get<2>(a), std::get<2>(b));
      }
    };
  }

  // Check the type used for ids
  template <class Grid, class IdSet>
  void checkIdType(const Grid& grid, const IdSet& idSet)
  {
    auto gridView = grid.leafGridView();
    auto begin = gridView.template begin<0>();

    // skip test for empty grids
    if( begin == gridView.template end<0>() )
      return ;

    // Get some entity for testing
    auto entity = *begin;

    // The IdSet class exports a type IdSet, and this is the type actually
    // used for return values of the id class
    static_assert(std::is_same<typename IdSet::IdType, decltype(idSet.id(entity))>::value,
                  "IdSet::IdType does not match the return value of idSet.id()");
    static_assert(std::is_same<typename IdSet::IdType, decltype(idSet.subId(entity,0,0))>::value,
                  "IdSet::IdType does not match the return value of idSet.subId()");

    // Get some id for testing
    auto id = idSet.id(entity);

    // Test for operator<
    if (id<id)
      DUNE_THROW(GridError, "operator< does not implement an ordering");

    // Test for operator==
    if (!(id==id))
      DUNE_THROW(GridError, "operator== is not symmetric");

    // Test for operator!=
    if (id!=id)
      DUNE_THROW(GridError, "operator!= is not the negation of operator==");

    // Test for the default constructor
    typename IdSet::IdType someId;

    // Test for the copy constructor
    typename IdSet::IdType someOtherId(idSet.id(entity));

    // copy assignment
    someId = someOtherId;

    // Output operator
    std::cout << someId << std::endl;
  }

  /** \brief Check that ids identify entities
   *
   * Two entities with the same id must have the same center, up to the
   * periodic map. Conversely, a periodic copy must have the same id as all
   * other entities on the same level whose centers are mapped to the same
   * point (this part is only active for a non-trivial periodic map).
   */
  template <class Grid, class IdSet>
  void checkInjectivity(const Grid& grid, const IdSet& idSet,
                        const GridCheck::PeriodicMap<Grid>& periodicMap = GridCheck::noPeriodicMap<Grid>())
  {
    const int dim = Grid::dimension;

    using IdType = typename IdSet::IdType;
    using GlobalCoordinate = typename Grid::template Codim<0>::Geometry::GlobalCoordinate;

    std::map<IdType,GlobalCoordinate> idContainer;

    const auto tolerance = GridCheck::coordinateTolerance(grid.leafGridView());

    // Loop over all entities of all codimensions on all grid levels
    for (int i=0; i<=grid.maxLevel(); ++i)
    {
      const auto gridView = grid.levelGridView(i);

      // ids of the entities on this level, by codim, geometry type, and mapped center
      // the flag records whether the id belongs to a periodic copy
      std::map<std::tuple<int,GeometryType,GlobalCoordinate>, std::pair<IdType,bool>,
               GridCheck::Impl::PositionLess> positionContainer(GridCheck::Impl::PositionLess{ tolerance });

      for (const auto& element: elements(gridView))
      {
        // Loop over all faces of all codimensions
        Hybrid::forEach( std::make_index_sequence< dim+1 >{},
                         [ & ]( auto codim )
        {
          if constexpr( Capabilities::hasEntity<Grid, codim>::v )
          {
            for (size_t face=0; face<element.subEntities(codim); ++face)
            {
              const auto& entity = element.template subEntity<codim>(face);

              auto id = idSet.id(entity);
              const auto center = entity.geometry().center();
              const auto mappedCenter = periodicMap(center);

              // Has the same id already been used by a different entity?
              if (idContainer.find(id) != idContainer.end())
              {
                // Yes.  Then either we have seen the same entity before, or we are now
                // on the copy of an entity we have seen before.  In either case we must
                // have the same entity center, up to the periodic map.
                // CAVEAT: This last reasoning does not hold if the grid uses parametrized
                // elements or parametrized boundaries.
                if (! GridCheck::coordinatesEqual(mappedCenter, idContainer[id], tolerance))
                  DUNE_THROW(GridError, "IdSet is not injective: entities at "
                             << idContainer[id] << " and " << mappedCenter
                             << " (codim " << codim << ") have the same id " << id);
              }
              else
              {
                idContainer[id] = mappedCenter;
              }

              // Is this a periodic copy of an entity on this level and process (or vice versa)?
              // Then both must have the same id. This only sees local coincidences; copies of
              // entities stored on other processes are checked by check_communication_correctness.
              const bool isCopy = GridCheck::isPeriodicCopy(periodicMap, center, tolerance);
              const auto key = std::make_tuple(int(codim), entity.type(), mappedCenter);
              auto it = positionContainer.find(key);
              if (it == positionContainer.end())
                positionContainer.emplace(key, std::make_pair(id, isCopy));
              else
              {
                if ((isCopy || it->second.second) && it->second.first != id)
                  DUNE_THROW(GridError, "Periodic copies have different ids: entity at "
                             << center << " (codim " << codim << ") has id " << id
                             << ", another copy of it has id " << it->second.first);
                it->second.second = it->second.second || isCopy;
              }

              // While we are here: Do subEntity.id and subId return the same value?
              if (id != idSet.subId(element,face,codim))
                DUNE_THROW(GridError, "subEntity.id and subId do not return the same value!");
            }
          }
        });
      }
    }
  }


  /** \brief Run all available tests for a given IdSet
   *
   * \param periodicMap maps each coordinate to its representative in the periodic
   *                    unit cell, see dune/grid/test/periodicmap.hh
   */
  template <class Grid, class IdSet>
  void checkIdSet ( const Grid &grid, const IdSet& idSet,
                    const GridCheck::PeriodicMap<Grid>& periodicMap = GridCheck::noPeriodicMap<Grid>())
  {
    checkIdType(grid,idSet);
    checkInjectivity(grid, idSet, periodicMap);
  }

} // end namespace Dune

#endif // DUNE_GRID_TEST_CHECKIDSET_HH
