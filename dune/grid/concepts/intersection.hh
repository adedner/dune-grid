// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_INTERSECTION_HH
#define DUNE_GRID_CONCEPTS_INTERSECTION_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/geometry.hh>

namespace Dune::Concept {

/**
 * @brief Model of an intersection
 * @ingroup GridConcepts
 * @details Dune::Grid::Intersection is a template for this model
 */
template<class I>
concept Intersection = requires(const I& i, typename I::LocalCoordinate local)
{
  requires EntityGeneral<typename I::Entity>;
  requires Geometry<typename I::Geometry>;
  requires Geometry<typename I::LocalGeometry>;
  typename I::ctype;
  { I::mydimension                  } -> std::convertible_to<int                            >;
  { I::dimensionworld               } -> std::convertible_to<int                            >;
  { i.boundary()                    } -> std::convertible_to<bool                           >;
  { i.boundarySegmentIndex()        } -> std::convertible_to<size_t                         >;
  { i.neighbor()                    } -> std::convertible_to<bool                           >;
  { i.inside()                      } -> std::convertible_to<typename I::Entity             >;
  { i.outside()                     } -> std::convertible_to<typename I::Entity             >;
  { i.conforming()                  } -> std::convertible_to<bool                           >;
  { i.geometryInInside()            } -> std::convertible_to<typename I::LocalGeometry      >;
  { i.geometryInOutside()           } -> std::convertible_to<typename I::LocalGeometry      >;
  { i.geometry()                    } -> std::convertible_to<typename I::Geometry           >;
  { i.type()                        } -> std::convertible_to<Dune::GeometryType             >;
  { i.indexInInside()               } -> std::convertible_to<int                            >;
  { i.indexInOutside()              } -> std::convertible_to<int                            >;
  { i.outerNormal(local)            } -> std::convertible_to<typename I::GlobalCoordinate   >;
  { i.integrationOuterNormal(local) } -> std::convertible_to<typename I::GlobalCoordinate   >;
  { i.unitOuterNormal(local)        } -> std::convertible_to<typename I::GlobalCoordinate   >;
  { i.centerUnitOuterNormal()       } -> std::convertible_to<typename I::GlobalCoordinate   >;
  { i==i                            } -> std::convertible_to<bool                           >;
  { i!=i                            } -> std::convertible_to<bool                           >;
  requires std::default_initializable<I>;
  requires std::copyable<I>;
};

} // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_INTERSECTION_HH
