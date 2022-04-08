// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRIDLEVELITERATOR_HH
#define DUNE_VIRTUALIZEDGRIDLEVELITERATOR_HH

#include <dune/grid/common/gridenums.hh>

/** \file
 * \brief The VirtualizedGridLevelIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup VirtualizedGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class VirtualizedGridLevelIterator
  {
  public:
    enum {codimension = codim};

    typedef typename GridImp::template Codim<codim>::Entity Entity;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual void increment () = 0;
      virtual Entity dereference () const = 0;
      virtual bool equals ( const VirtualizedGridLevelIterator& i ) const = 0;
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
        return VirtualizedGridEntity<codim, GridImp::dimension, GridImp> ( std::move(*impl()) );
      }

      virtual bool equals( const VirtualizedGridLevelIterator& i ) const override
      {
        return impl() == dynamic_cast<Implementation<I>&>(*i.impl_).impl();
      }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I impl_;
    };
    // VIRTUALIZATION END


  public:
    template< class ImplLevelIterator >
    explicit VirtualizedGridLevelIterator(ImplLevelIterator&& implLevelIterator)
    : impl_( new Implementation<ImplLevelIterator>( std::forward<ImplLevelIterator>(implLevelIterator) ) )
    {}

    VirtualizedGridLevelIterator(const VirtualizedGridLevelIterator& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLevelIterator ( VirtualizedGridLevelIterator && ) = default;

    VirtualizedGridLevelIterator& operator=(const VirtualizedGridLevelIterator& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    //! prefix increment
    void increment() {
      impl_->increment();
    }

    //! dereferencing
    Entity dereference() const {
      return impl_->dereference();
    }

    //! equality
    bool equals(const VirtualizedGridLevelIterator& i) const {
      return impl_->equals(i);
    }

  private:
    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune

#endif
