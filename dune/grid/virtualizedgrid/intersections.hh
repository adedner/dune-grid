// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRID_INTERSECTIONS_HH
#define DUNE_VIRTUALIZEDGRID_INTERSECTIONS_HH

#include "virtualizedgridleafiterator.hh"
#include <dune/grid/virtualizedgrid/virtualizedgrvirtualized.hh>

/** \file
 * \brief The VirtualizedGridLeafIntersection and VirtualizedGridLevelIntersection classes
 */

namespace Dune {


  /** \brief An intersection with a leaf neighbor element
   * \ingroup VirtualizedGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class VirtualizedGridLeafIntersection
  {

    friend class VirtualizedGridLeafIntersectionIterator<GridImp>;

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals(const VirtualizedGridLeafIntersection& other) const = 0;
      virtual Entity inside() const = 0;
      virtual Entity outside() const = 0;
      virtual bool boundary () const = 0;
      virtual NormalVector centerUnitOuterNormal () const = 0;
      virtual bool neighbor () const = 0;
      virtual size_t boundarySegmentIndex() const = 0;
      virtual bool conforming () const = 0;
      virtual GeometryType type () const = 0;
      virtual LocalGeometry geometryInInside () const = 0;
      virtual LocalGeometry geometryInOutside () const = 0;
      virtual Geometry geometry () const = 0;
      virtual int indexInInside () const = 0;
      virtual int indexInOutside () const = 0;
      virtual FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const = 0;
      virtual FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const = 0;
      virtual FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool equals(const VirtualizedGridLeafIntersection& other) const override { return impl().equals(other); }
      virtual Entity inside() const override { return impl().inside(); }
      virtual Entity outside() const override { return impl().outside(); }
      virtual bool boundary () const override { return impl().boundary(); }
      virtual NormalVector centerUnitOuterNormal () const override { return impl().centerUnitOuterNormal(); }
      virtual bool neighbor () const override { return impl().neighbor(); }
      virtual size_t boundarySegmentIndex() const override { return impl().boundarySegmentIndex(); }
      virtual bool conforming () const override { return impl().conforming(); }
      virtual GeometryType type () const override { return impl().type(); }
      virtual LocalGeometry geometryInInside () const override { return impl().geometryInInside(); }
      virtual LocalGeometry geometryInOutside () const override { return impl().geometryInOutside(); }
      virtual Geometry geometry () const override { return impl().geometry(); }
      virtual int indexInInside () const override { return impl().indexInInside(); }
      virtual int indexInOutside () const override { return impl().indexInOutside(); }
      virtual FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const override { return impl().outerNormal(local); }
      virtual FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const override { return impl().integrationOuterNormal(local); }
      virtual FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const override { return impl().unitOuterNormal(local); }

    private:
      const auto &impl () const { return impl_; }
      const I& impl_;
    };
    // VIRTUALIZATION END

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    VirtualizedGridLeafIntersection()
    {}

    template< class ImplIntersection >
    VirtualizedGridLeafIntersection(const GridImp* virtualizedGrid,
                                    const ImplIntersection& implIntersection)
      : virtualizedGrid_(virtualizedGrid),
        impl_( new Implementation(implIntersection) )
    {}

    bool equals(const VirtualizedGridLeafIntersection& other) const
    {
      return impl_ == other.impl_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return VirtualizedGridEntity<0,dim,GridImp>(virtualizedGrid_, impl_->inside());
    }


    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {
      return VirtualizedGridEntity<0,dim,GridImp>(virtualizedGrid_, impl_->outside());
    }


    //! return true if intersection is with boundary.
    bool boundary () const {
      return impl_->boundary();
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return impl_->centerUnitOuterNormal();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return impl_->neighbor();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      return impl_->boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      return impl_->conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return impl_->type();
    }


    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( impl_->geometryInInside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( impl_->geometryInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( impl_->geometry() );
    }


    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return impl_->indexInInside();
    }


    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      return impl_->indexInOutside();
    }


    //! return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return impl_->outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return impl_->integrationOuterNormal(local);
    }

    //! return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return impl_->unitOuterNormal(local);
    }

  private:
    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;
  };




  //! \todo Please doc me !
  template<class GridImp>
  class VirtualizedGridLevelIntersection
  {

    friend class VirtualizedGridLevelIntersectionIterator<GridImp>;

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals(const VirtualizedGridLeafIntersection& other) const = 0;
      virtual Entity inside() const = 0;
      virtual Entity outside() const = 0;
      virtual bool boundary () const = 0;
      virtual NormalVector centerUnitOuterNormal () const = 0;
      virtual bool neighbor () const = 0;
      virtual size_t boundarySegmentIndex() const = 0;
      virtual bool conforming () const = 0;
      virtual GeometryType type () const = 0;
      virtual LocalGeometry geometryInInside () const = 0;
      virtual LocalGeometry geometryInOutside () const = 0;
      virtual Geometry geometry () const = 0;
      virtual int indexInInside () const = 0;
      virtual int indexInOutside () const = 0;
      virtual FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const = 0;
      virtual FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const = 0;
      virtual FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool equals(const VirtualizedGridLeafIntersection& other) const override { return impl().equals(other); }
      virtual Entity inside() const override { return impl().inside(); }
      virtual Entity outside() const override { return impl().outside(); }
      virtual bool boundary () const override { return impl().boundary(); }
      virtual NormalVector centerUnitOuterNormal () const override { return impl().centerUnitOuterNormal(); }
      virtual bool neighbor () const override { return impl().neighbor(); }
      virtual size_t boundarySegmentIndex() const override { return impl().boundarySegmentIndex(); }
      virtual bool conforming () const override { return impl().conforming(); }
      virtual GeometryType type () const override { return impl().type(); }
      virtual LocalGeometry geometryInInside () const override { return impl().geometryInInside(); }
      virtual LocalGeometry geometryInOutside () const override { return impl().geometryInOutside(); }
      virtual Geometry geometry () const override { return impl().geometry(); }
      virtual int indexInInside () const override { return impl().indexInInside(); }
      virtual int indexInOutside () const override { return impl().indexInOutside(); }
      virtual FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const override { return impl().outerNormal(local); }
      virtual FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const override { return impl().integrationOuterNormal(local); }
      virtual FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const override { return impl().unitOuterNormal(local); }

    private:
      const auto &impl () const { return impl_; }
      const I& impl_;
    };
    // VIRTUALIZATION END

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    VirtualizedGridLevelIntersection()
    {}

    template< class ImplIntersection >
    VirtualizedGridLevelIntersection(const GridImp* virtualizedGrid,
                                     const ImplIntersection& implIntersection)
      : virtualizedGrid_(virtualizedGrid),
        impl_( new Implementation(implIntersection) )
    {}

    bool equals(const VirtualizedGridLevelIntersection& other) const
    {
      return impl_ == other.impl_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return VirtualizedGridEntity<0,dim,GridImp>(virtualizedGrid_, impl_->inside());
    }


    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {
      return VirtualizedGridEntity<0,dim,GridImp>(virtualizedGrid_, impl_->outside());
    }


    /** \brief return true if intersection is with boundary.
     */
    bool boundary () const {
      return impl_->boundary();
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return impl_->centerUnitOuterNormal();
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return impl_->neighbor();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      return impl_->boundarySegmentIndex();
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      return impl_->conforming();
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return impl_->type();
    }


    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( impl_->geometryInInside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( impl_->geometryInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( impl_->geometry() );
    }


    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return impl_->indexInInside();
    }


    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      return impl_->indexInOutside();
    }


    //! return outer normal
    FieldVector<ctype, dimworld> outerNormal (const FieldVector<ctype, dim-1>& local) const {
      return impl_->outerNormal(local);
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, dimworld> integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const {
      return impl_->integrationOuterNormal(local);
    }

    //! return unit outer normal
    FieldVector<ctype, dimworld> unitOuterNormal (const FieldVector<ctype, dim-1>& local) const {
      return impl_->unitOuterNormal(local);
    }

  private:
    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune

#endif
