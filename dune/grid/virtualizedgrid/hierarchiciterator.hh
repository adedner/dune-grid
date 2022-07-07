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
  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual void increment () = 0;
      virtual Entity dereference () const = 0;
      virtual bool equals ( const VirtualizedGridHierarchicIterator<GridImp>& i ) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::forward<I>(i) ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual void increment() override { ++impl(); }

      virtual Entity dereference() const override
      {
        return VirtualizedGridEntity<codimension, GridImp::dimension, GridImp> ( *impl() );
      }

      virtual bool equals( const VirtualizedGridHierarchicIterator<GridImp>& i ) const override
      {
        return impl() == static_cast<Implementation<I>&>(*i.impl_).impl();
      }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I impl_;
    };
    // VIRTUALIZATION END

  public:
    template< class ImplHierarchicIterator >
    explicit VirtualizedGridHierarchicIterator(ImplHierarchicIterator&& implHierarchicIterator)
    : impl_( new Implementation<ImplHierarchicIterator>( std::forward<ImplHierarchicIterator>( implHierarchicIterator ) ) )
    {}

    VirtualizedGridHierarchicIterator(const VirtualizedGridHierarchicIterator& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridHierarchicIterator ( VirtualizedGridHierarchicIterator && ) = default;

    VirtualizedGridHierarchicIterator& operator=(const VirtualizedGridHierarchicIterator& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    void increment()
    {
      impl_->increment();
    }

    Entity dereference() const {
      return impl_->dereference();
    }

    bool equals(const VirtualizedGridHierarchicIterator& i) const {
      return impl_->equals(i);
    }

    std::unique_ptr<Interface> impl_;
  };


}  // end namespace Dune

#endif
