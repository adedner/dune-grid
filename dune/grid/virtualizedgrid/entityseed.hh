// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH
#define DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH

/**
 * \file
 * \brief The VirtualizedGridEntitySeed class
 */


namespace Dune {


  /**
   * \brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup VirtualizedGrid
   *
   */
  template<int codim, class GridImp>
  class VirtualizedGridEntitySeed
  {
  protected:

    // Entity type of the grid
    typedef typename GridImp::Traits::template Codim<codim>::Entity Entity;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool isValid () const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool isValid() const override { impl().increment(); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END

  public:

    enum {codimension = codim};

    /**
     * \brief Construct an empty (i.e. isValid() == false) seed.
     */
    VirtualizedGridEntitySeed()
    {}

    /**
     * \brief Create EntitySeed from hostgrid Entity
     *
     * We call hostEntity.seed() directly in the constructor
     * of VirtualizedGridEntitySeed to allow for return value optimization.
     */
    VirtualizedGridEntitySeed(const Entity& entity)
    : impl_( new Implementation( entity.impl().impl_ ) )
    {}

    /**
     * \brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid() const
    {
      return impl_->isValid();
    }
  private:

    std::unique_ptr<Interface> impl_;
  };

} // namespace Dune

#endif  // #define DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH
