// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRIDHIERITERATOR_HH
#define DUNE_VIRTUALIZEDGRIDHIERITERATOR_HH

/** \file
 * \brief The VirtualizedGridHierarchicIterator class
 */

namespace Dune {


  //**********************************************************************
  //
  /** \brief Iterator over the descendants of an entity.
   * \ingroup VirtualizedGrid
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HierarchicIterator,
     starting from a given entity.
   */
  template<class GridImp>
  class VirtualizedGridHierarchicIterator
  {

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual void increment () = 0;
      virtual Entity dereference () = 0;
      virtual bool equals ( const VirtualizedGridLeafIterator& i ) = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual void increment() override { impl().increment(); }
      virtual Entity dereference() const override { return impl().dereference(); }
      virtual bool equals( const VirtualizedGridHierarchicIterator& i ) const override { return impl().equals(i); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END

  public:

    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! the default Constructor
    explicit VirtualizedGridHierarchicIterator(const GridImp* virtualizedGrid, const Entity& startEntity, int maxLevel) :
      virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation( startEntity.impl().impl_.hbegin(maxLevel) )
    {}


    //! \todo Please doc me !
    explicit VirtualizedGridHierarchicIterator(const GridImp* virtualizedGrid, const Entity& startEntity, int maxLevel, [[maybe_unused]] bool endDummy) :
      virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation( startEntity.impl().impl_.hend(maxLevel) )
    {}


    //! \todo Please doc me !
    void increment()
    {
      impl_->increment();
    }

    //! dereferencing
    Entity dereference() const {
      return Entity{{virtualizedGrid_,*hostHierarchicIterator_}};
    }

    //! equality
    bool equals(const VirtualizedGridHierarchicIterator& i) const {
      return impl_ == i.impl_;
    }

  private:
    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;

  };


}  // end namespace Dune

#endif
