// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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
      virtual bool equals( const VirtualizedGridLeafIterator& i ) const override { return impl().equals(i); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END


  public:
    explicit VirtualizedGridLeafIterator(const GridImp* virtualizedGrid)
    : virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation( virtualizedGrid->impl_->leafGridView().template begin<codim,pitype>() ) )
    {}

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param virtualizedGrid  pointer to grid instance
     */
    explicit VirtualizedGridLeafIterator(const GridImp* virtualizedGrid, [[maybe_unused]] bool endDummy)
    : virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation( virtualizedGrid->impl_->leafGridView().template end<codim,pitype>() ) )
    {}


    //! prefix increment
    void increment() {
      impl_->increment();
    }

    //! dereferencing
    Entity dereference() const {
      return Entity{{virtualizedGrid_, impl_->dereference()}};
    }

    //! equality
    bool equals(const VirtualizedGridLeafIterator& i) const {
      return impl_ == i.impl_;
    }

  private:
    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;

  };


}  // namespace Dune

#endif
