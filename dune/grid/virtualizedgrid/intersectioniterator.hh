// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRID_INTERSECTIONITERATOR_HH
#define DUNE_VIRTUALIZEDGRID_INTERSECTIONITERATOR_HH

#include "virtualizedgridintersections.hh"
#include "virtualizedgrvirtualized.hh"

#include <dune/grid/common/intersection.hh>

/** \file
 * \brief The VirtualizedGridLeafIntersectionIterator and VirtualizedGridLevelIntersectionIterator classes
 */

namespace Dune {

  /** \brief Iterator over all element neighbors
   * \ingroup VirtualizedGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class VirtualizedGridLeafIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

  public:
    typedef Dune::Intersection<const GridImp, Dune::VirtualizedGridLeafIntersection<GridImp> > Intersection;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals ( const VirtualizedGridLeafIterator& i ) = 0;
      virtual void increment () = 0;
      virtual Intersection dereference () = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool equals( const VirtualizedGridLeafIntersectionIterator& i ) const override { return impl().equals(i); }
      virtual void increment() override { impl().increment(); }
      virtual Intersection dereference() const override { return impl().dereference(); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END


  public:
    VirtualizedGridLeafIntersectionIterator()
    {}

    template< class ImplLeafIntersectionIterator >
    VirtualizedGridLeafIntersectionIterator(const GridImp* virtualizedGrid,
                                            const ImplLeafIntersectionIterator& implIterator)
      : virtualizedGrid_(virtualizedGrid),
        impl_( new Implementation(implIterator) )
    {}

    //! equality
    bool equals(const VirtualizedGridLeafIntersectionIterator& other) const {
      return impl_ == other.impl_;
    }


    //! prefix increment
    void increment() {
      impl_->increment();
    }

    //! \brief dereferencing
    Intersection dereference() const {
      return VirtualizedGridLeafIntersection<GridImp>(virtualizedGrid_, impl_->dereference());
    }

  private:
    //**********************************************************
    //  private data
    //**********************************************************

    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;
  };




  //! \todo Please doc me !
  template<class GridImp>
  class VirtualizedGridLevelIntersectionIterator
  {
    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

  public:

    typedef Dune::Intersection<const GridImp, Dune::VirtualizedGridLevelIntersection<GridImp> > Intersection;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals ( const VirtualizedGridLevelIntersectionIterator<GridImp>& i ) = 0;
      virtual void increment () = 0;
      virtual Intersection dereference () = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool equals( const VirtualizedGridLevelIntersectionIterator<GridImp>& i ) const override { return impl().equals(i); }
      virtual void increment() override { impl().increment(); }
      virtual Intersection dereference() const override { return impl().dereference(); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END

  public:
    VirtualizedGridLevelIntersectionIterator()
    {}

    template< class ImplLevelIntersectionIterator >
    VirtualizedGridLevelIntersectionIterator(const GridImp* virtualizedGrid,
                                             const ImplLevelIntersectionIterator& implIterator)
      : virtualizedGrid_(virtualizedGrid),
        impl_( new Implementation(implIterator) )
    {}

    //! equality
    bool equals(const VirtualizedGridLevelIntersectionIterator<GridImp>& other) const {
      return impl_ == other.impl_;
    }


    //! prefix increment
    void increment() {
      impl_->increment();
    }

    //! \brief dereferencing
    Intersection dereference() const {
      return VirtualizedGridLevelIntersection<GridImp>(virtualizedGrid_, impl_->dereference());
    }

  private:
    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune

#endif
