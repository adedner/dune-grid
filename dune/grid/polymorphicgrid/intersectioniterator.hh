// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRID_INTERSECTIONITERATOR_HH
#define DUNE_VIRTUALIZEDGRID_INTERSECTIONITERATOR_HH

#include "intersections.hh"
#include "entity.hh"

#include <dune/grid/common/intersection.hh>
#include <dune/grid/polymorphicgrid/common/iterator.hh>

/** \file
 * \brief The PolymorphicIntersectionIterator class
 */

namespace Dune {

  /** \brief Iterator over all element neighbors
   * \ingroup Polymorphic
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class PolymorphicIntersectionIterator
      : public Polymorphic::IteratorDefinition<PolymorphicIntersectionIterator<GridImp>,PolymorphicIntersection<GridImp>>::Base
  {
    using Self = PolymorphicIntersectionIterator;
    using IntersectionImp = PolymorphicIntersection<const GridImp>;
    using Definition = Polymorphic::IteratorDefinition<Self,IntersectionImp>;
    using Base = typename Definition::Base;

  public:

    enum {dim=GridImp::dimension};
    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    using ctype = typename GridImp::ctype;

  public:
    using Intersection = Dune::Intersection<const GridImp, IntersectionImp>;

    PolymorphicIntersectionIterator () = default;

    template <class Impl, disableCopyMove<PolymorphicIntersectionIterator,Impl> = 0>
    PolymorphicIntersectionIterator (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    //! equality
    bool equals (const Self& other) const
    {
      return this->asInterface().equals(other);
    }

    //! prefix increment
    void increment ()
    {
      this->asInterface().increment();
    }

    //! dereferencing
    Intersection dereference () const
    {
      return Intersection{ this->asInterface().dereference() };
    }
  };

}  // namespace Dune

#endif
