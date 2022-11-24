// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRIDLEAFITERATOR_HH
#define DUNE_VIRTUALIZEDGRIDLEAFITERATOR_HH

#include <dune/grid/common/gridenums.hh>

/** \file
 * \brief The VirtualizedGridLeafIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup VirtualizedGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class VirtualizedGridLeafIterator
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
      virtual bool equals ( const VirtualizedGridLeafIterator& i ) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::forward<I>(i) ) {}
      Implementation *clone() const override { return new Implementation( *this ); }

      void increment() override { ++impl(); }

      Entity dereference() const override
      {
        return VirtualizedGridEntity<codim, GridImp::dimension, GridImp> ( std::move( *impl() ) );
      }

      bool equals( const VirtualizedGridLeafIterator& i ) const override
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
    VirtualizedGridLeafIterator () = default;

    template< class ImplLeafIterator >
    explicit VirtualizedGridLeafIterator(ImplLeafIterator&& implLeafIterator)
    : impl_( new Implementation<ImplLeafIterator>( std::forward<ImplLeafIterator>(implLeafIterator) ) )
    {}

    VirtualizedGridLeafIterator(const VirtualizedGridLeafIterator& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLeafIterator ( VirtualizedGridLeafIterator && ) = default;

    VirtualizedGridLeafIterator& operator=(const VirtualizedGridLeafIterator& other)
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
    bool equals(const VirtualizedGridLeafIterator& i) const {
      return impl_->equals(i);
    }

    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune

#endif
