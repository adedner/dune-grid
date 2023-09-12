// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRID_INTERSECTIONS_HH
#define DUNE_VIRTUALIZEDGRID_INTERSECTIONS_HH

#include "entity.hh"

/** \file
 * \brief The VirtualizedGridIntersection class
 */

namespace Dune {

  template<class GridImp, class Derived>
  struct VirtualizedGridIntersectionDefinition
  {
    enum {dim=GridImp::dimension};
    enum {dimworld=GridImp::dimensionworld};

    using Geometry = typename GridImp::template Codim<1>::Geometry;
    using GeometryImp = VirtualizedGridGeometry<dim-1, Geometry::coorddimension, GridImp>;
    using LocalGeometry = typename GridImp::template Codim<1>::LocalGeometry;
    using LocalGeometryImp = VirtualizedGridGeometry<dim-1, dim, GridImp>;
    using Entity = typename GridImp::template Codim<0>::Entity;
    using EntityImp = VirtualizedGridEntity<0, dim, GridImp>;
    using LocalCoordinate = typename Geometry::LocalCoordinate;
    using NormalVector = FieldVector<typename GridImp::ctype, dimworld>;

    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool equals (const Derived& other) const = 0;
      virtual Entity inside () const = 0;
      virtual Entity outside () const = 0;
      virtual bool boundary () const = 0;
      virtual bool neighbor () const = 0;
      virtual size_t boundarySegmentIndex () const = 0;
      virtual bool conforming () const = 0;
      virtual GeometryType type () const = 0;
      virtual LocalGeometry geometryInInside () const = 0;
      virtual LocalGeometry geometryInOutside () const = 0;
      virtual Geometry geometry () const = 0;
      virtual int indexInInside () const = 0;
      virtual int indexInOutside () const = 0;
      virtual NormalVector outerNormal (const LocalCoordinate& local) const = 0;
      virtual NormalVector integrationOuterNormal (const LocalCoordinate& local) const = 0;
      virtual NormalVector unitOuterNormal (const LocalCoordinate& local) const = 0;
      virtual NormalVector centerUnitOuterNormal () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      using Wrapped = typename Wrapper::Wrapped;

      bool equals (const Derived& other) const final { return this->get()() == Polymorphic::asWrapped<Wrapped>(other); }
      Entity inside () const final { return EntityImp{this->get().inside()}; }
      Entity outside () const final { return EntityImp{this->get().outside()}; }
      bool boundary () const final { return this->get().boundary(); }
      bool neighbor () const final { return this->get().neighbor(); }
      size_t boundarySegmentIndex () const final { return this->get().boundarySegmentIndex(); }
      bool conforming () const final { return this->get().conforming(); }
      GeometryType type () const final { return this->get().type(); }
      LocalGeometry geometryInInside () const final { return LocalGeometryImp{this->get().geometryInInside()}; }
      LocalGeometry geometryInOutside () const final { return LocalGeometryImp{this->get().geometryInOutside()}; }
      Geometry geometry () const final { return GeometryImp{this->get().geometry()}; }
      int indexInInside () const final { return this->get().indexInInside(); }
      int indexInOutside () const final { return this->get().indexInOutside(); }
      NormalVector outerNormal (const LocalCoordinate& local) const final { return this->get().outerNormal(local); }
      NormalVector integrationOuterNormal (const LocalCoordinate& local) const final { return this->get().integrationOuterNormal(local); }
      NormalVector unitOuterNormal (const LocalCoordinate& local) const override { return this->get().unitOuterNormal(local); }
      NormalVector centerUnitOuterNormal () const final { return this->get().centerUnitOuterNormal(); }
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };

  /** \brief An intersection with a leaf neighbor element
   * \ingroup VirtualizedGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class VirtualizedGridIntersection
  {
    using Self = VirtualizedGridIntersection;
    using Definition = VirtualizedGridIntersectionDefinition<GridImp,Self>;
    using Base = typename Definition::Base;

    friend class VirtualizedGridIntersectionIterator<GridImp>;

  public:
    enum {dim=GridImp::dimension};
    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    using ctype = typename GridImp::ctype;

  public:
    using Geometry = typename GridImp::template Codim<1>::Geometry;
    using LocalGeometry = typename GridImp::template Codim<1>::LocalGeometry;
    using Entity = typename GridImp::template Codim<0>::Entity;
    using LocalCoordinate = typename Geometry::LocalCoordinate;
    using NormalVector = FieldVector<ctype, dimworld>;

  public:
    VirtualizedGridIntersection () = default;

    template <class Impl, disableCopyMove<VirtualizedGridIntersection,Impl> = 0>
    VirtualizedGridIntersection (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    bool equals (const Self& other) const
    {
      return this->asInterface().equals(other);
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside () const
    {
      return this->asInterface().inside();
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside () const
    {
      return this->asInterface().outside();
    }

    //! return true if intersection is with boundary.
    bool boundary () const
    {
      return this->asInterface().boundary();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const
    {
      return this->asInterface().neighbor();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const
    {
      return this->asInterface().boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const
    {
      return this->asInterface().conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const
    {
      return this->asInterface().type();
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return this->asInterface().geometryInInside();
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return this->asInterface().geometryInOutside();
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return this->asInterface().geometry();
    }

    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return this->asInterface().indexInInside();
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const
    {
      return this->asInterface().indexInOutside();
    }

    //! return outer normal
    NormalVector outerNormal (const LocalCoordinate& local) const
    {
      return this->asInterface().outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    NormalVector integrationOuterNormal (const LocalCoordinate& local) const
    {
      return this->asInterface().integrationOuterNormal(local);
    }

    //! return unit outer normal
    NormalVector unitOuterNormal (const LocalCoordinate& local) const
    {
      return this->asInterface().unitOuterNormal(local);
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const
    {
      return this->asInterface().centerUnitOuterNormal();
    }
  };

}  // namespace Dune

#endif
