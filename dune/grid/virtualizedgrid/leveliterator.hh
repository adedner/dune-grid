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
   * \ingroup VirtualizedGrid
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
      virtual Entity dereference () = 0;
      virtual bool equals ( const VirtualizedGridLevelIterator& i ) = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual void increment() override { impl().increment(); }
      virtual Entity dereference() const override { return impl().dereference(); }
      virtual bool equals( const VirtualizedGridLevelIterator& i ) const override { return impl().equals(i); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END


  public:
    explicit VirtualizedGridLevelIterator(const GridImp* virtualizedGrid, int level)
    : virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation( virtualizedGrid->impl_->levelGridView(level).template begin<codim,pitype>() ) )
    {}


    /** \brief Constructor which create the end iterator
        \param endDummy      Here only to distinguish it from the other constructor
        \param virtualizedGrid  pointer to VirtualizedGrid instance
        \param level         grid level on which the iterator shall be created
     */
    explicit VirtualizedGridLevelIterator(const GridImp* virtualizedGrid, int level, [[maybe_unused]] bool endDummy)
    : virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation( virtualizedGrid->impl_->levelGridView(level).template end<codim,pitype>() ) )
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
    bool equals(const VirtualizedGridLevelIterator& i) const {
      return impl_ == i.impl_;
    }

  private:
    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune

#endif
