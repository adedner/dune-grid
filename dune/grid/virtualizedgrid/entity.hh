// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRVIRTUALIZED_HH
#define DUNE_VIRTUALIZEDGRVIRTUALIZED_HH

/** \file
 * \brief The VirtualizedGridEntity class
 */

#include <dune/grid/common/grid.hh>

namespace Dune {


  // Forward declarations

  template<int codim, int dim, class GridImp>
  class VirtualizedGridEntity;

  template<int codim, class GridImp>
  class VirtualizedGridEntitySeed;

  template<int codim, PartitionIteratorType pitype, class GridImp>
  class VirtualizedGridLevelIterator;

  template<class GridImp>
  class VirtualizedGridLevelIntersectionIterator;

  template<class GridImp>
  class VirtualizedGridLeafIntersectionIterator;

  template<class GridImp>
  class VirtualizedGridHierarchicIterator;

  //**********************************************************************
  //
  // --VirtualizedGridEntity
  // --Entity
  //
  /** \brief The implementation of entities in a VirtualizedGrid
   *   \ingroup VirtualizedGrid
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int codim, int dim, class GridImp>
  class VirtualizedGridEntity :
    public EntityDefaultImplementation <codim, dim, GridImp, VirtualizedGridEntity>
  {

    template <class GridImp_>
    friend class VirtualizedGridLevelIndexSet;

    template <class GridImp_>
    friend class VirtualizedGridLeafIndexSet;

    template <class GridImp_>
    friend class VirtualizedGridLocalIdSet;

    template <class GridImp_>
    friend class VirtualizedGridGlobalIdSet;

  public:
    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;

  private:

    typedef typename GridImp::ctype ctype;

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals(const VirtualizedGridEntity<codim, dim, GridImp>& other) const = 0;
      virtual EntitySeed seed () const = 0;
      virtual int level () const = 0;
      virtual PartitionType partitionType () const = 0;
      virtual unsigned int subEntities (unsigned int cc) const = 0;
      virtual Geometry geometry () const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::forward<I>(i) ) {}

      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual bool equals(const VirtualizedGridEntity<codim, dim, GridImp>& other) const override
      {
        return impl() == static_cast<Implementation<I>&>(*other.impl_).impl();
      }

      virtual EntitySeed seed () const override
      {
        return VirtualizedGridEntitySeed<codim, GridImp>( impl().seed() );
      }

      virtual int level () const override { return impl().level(); }

      virtual PartitionType partitionType () const override { return impl().partitionType(); }

      virtual unsigned int subEntities (unsigned int cc) const override { return impl().subEntities(cc); }

      virtual Geometry geometry () const override
      {
        return Geometry( VirtualizedGridGeometry<dim-codim, dim, GridImp>( impl().geometry() ) );
      }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };
    // VIRTUALIZATION END

  public:

    VirtualizedGridEntity()
    {}

    template< class ImplGridEntity >
    VirtualizedGridEntity(ImplGridEntity&& implEntity)
    : impl_( new Implementation<const ImplGridEntity>( std::forward<ImplGridEntity>( implEntity ) ) )
    {}

    VirtualizedGridEntity(const VirtualizedGridEntity& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridEntity ( VirtualizedGridEntity && ) = default;

    VirtualizedGridEntity& operator=(const VirtualizedGridEntity& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    VirtualizedGridEntity& operator=( VirtualizedGridEntity&& ) = default;

    bool equals(const VirtualizedGridEntity& other) const
    {
      return impl_->equals(other);
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return impl_->seed();
    }

    //! level of this element
    int level () const {
      return impl_->level();
    }

    /** \brief The partition type for parallel computing
     */
    PartitionType partitionType () const {
      return impl_->partitionType();
    }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cc) const
    {
      return impl_->subEntities(cc);
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      return Geometry( impl_->geometry() );
    }

    std::unique_ptr<Interface> impl_;
  };




  //***********************
  //
  //  --VirtualizedGridEntity
  //
  //***********************
  /** \brief Specialization for codim-0-entities.
   * \ingroup VirtualizedGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   */
  template<int dim, class GridImp>
  class VirtualizedGridEntity<0, dim, GridImp> :
    public EntityDefaultImplementation<0, dim, GridImp, VirtualizedGridEntity>
  {
  public:
    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on this level
    typedef VirtualizedGridLevelIntersectionIterator<const GridImp> LevelIntersectionIterator;

    //! The Iterator over intersections on the leaf level
    typedef VirtualizedGridLeafIntersectionIterator<const GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef VirtualizedGridHierarchicIterator<const GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

  private:

    typedef typename GridImp::ctype ctype;

    // VIRTUALIZATION BEGIN
  public:
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals(const VirtualizedGridEntity<0, dim, GridImp>& other) const = 0;
      virtual bool hasFather () const = 0;
      virtual EntitySeed seed () const = 0;
      virtual int level () const = 0;
      virtual PartitionType partitionType () const = 0;
      virtual Geometry geometry () const = 0;
      virtual unsigned int subEntities (unsigned int cc) const = 0;
      virtual typename GridImp::template Codim<0>::Entity subEntity0 (int i) const = 0;
      virtual typename GridImp::template Codim<1>::Entity subEntity1 (int i) const = 0;
      virtual typename GridImp::template Codim<dim-1>::Entity subEntityDimMinus1 (int i) const = 0;
      virtual typename GridImp::template Codim<dim>::Entity subEntityDim (int i) const = 0;
      virtual LevelIntersectionIterator ilevelbegin () const = 0;
      virtual LevelIntersectionIterator ilevelend () const = 0;
      virtual LeafIntersectionIterator ileafbegin () const = 0;
      virtual LeafIntersectionIterator ileafend () const = 0;
      virtual bool isLeaf() const = 0;
      virtual typename GridImp::template Codim<0>::Entity father () const = 0;
      virtual LocalGeometry geometryInFather () const = 0;
      virtual HierarchicIterator hbegin (int maxLevel) const = 0;
      virtual HierarchicIterator hend (int maxLevel) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::forward<I>(i) ) {}

      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual bool equals(const VirtualizedGridEntity<0, dim, GridImp>& other) const override {
        return impl() == static_cast<Implementation<I>&>(*other.impl_).impl();
      }

      virtual bool hasFather () const override { return impl().hasFather(); }

      virtual EntitySeed seed () const override
      {
        return VirtualizedGridEntitySeed<0, GridImp>( impl().seed() );
      }

      virtual int level () const override { return impl().level(); }

      virtual PartitionType partitionType () const override { return impl().partitionType(); }

      virtual Geometry geometry () const override
      {
        return Geometry( VirtualizedGridGeometry<dim, dim, GridImp>( impl().geometry() ) );
      }

      virtual unsigned int subEntities (unsigned int cc) const override { return impl().subEntities(cc); }

      virtual typename GridImp::template Codim<0>::Entity subEntity0 (int i) const override
      {
        return VirtualizedGridEntity<0, dim, GridImp>( impl().template subEntity<0>(i) );
      }

      virtual typename GridImp::template Codim<1>::Entity subEntity1 (int i) const override
      {
        return VirtualizedGridEntity<1, dim, GridImp>( impl().template subEntity<1>(i) );
      }

      virtual typename GridImp::template Codim<dim-1>::Entity subEntityDimMinus1 (int i) const override
      {
        return VirtualizedGridEntity<dim-1, dim, GridImp>( impl().template subEntity<dim-1>(i) );
      }

      virtual typename GridImp::template Codim<dim>::Entity subEntityDim (int i) const override
      {
        return VirtualizedGridEntity<dim, dim, GridImp>( impl().template subEntity<dim>(i) );
      }

      virtual LevelIntersectionIterator ilevelbegin () const override
      {
        return VirtualizedGridLevelIntersectionIterator<const GridImp>( impl().impl().ilevelbegin() );
      }

      virtual LevelIntersectionIterator ilevelend () const override
      {
        return VirtualizedGridLevelIntersectionIterator<const GridImp>( impl().impl().ilevelend() );
      }

      virtual LeafIntersectionIterator ileafbegin () const override
      {
        return VirtualizedGridLeafIntersectionIterator<const GridImp>( impl().impl().ileafbegin() );
      }

      virtual LeafIntersectionIterator ileafend () const override
      {
        return VirtualizedGridLeafIntersectionIterator<const GridImp>( impl().impl().ileafend() );
      }

      virtual bool isLeaf() const override { return impl().isLeaf(); }

      virtual typename GridImp::template Codim<0>::Entity father () const override
      {
        return VirtualizedGridEntity<0, dim, GridImp>(impl().father());
      }

      virtual LocalGeometry geometryInFather () const override
      {
        return LocalGeometry( VirtualizedGridGeometry<dim, dim, GridImp>( impl().geometryInFather() ) );
      }

      virtual HierarchicIterator hbegin (int maxLevel) const override
      {
        return VirtualizedGridHierarchicIterator<const GridImp>(impl().hbegin(maxLevel));
      }

      virtual HierarchicIterator hend (int maxLevel) const override
      {
        return VirtualizedGridHierarchicIterator<const GridImp>(impl().hend(maxLevel));
      }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };
    // VIRTUALIZATION END

    VirtualizedGridEntity()
    {}

    template< class ImplGridEntity >
    VirtualizedGridEntity(ImplGridEntity&& implEntity)
    : impl_( new Implementation<const ImplGridEntity>( std::forward<ImplGridEntity>( implEntity ) ) )
    {}

    VirtualizedGridEntity(const VirtualizedGridEntity& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridEntity ( VirtualizedGridEntity && ) = default;

    VirtualizedGridEntity& operator=(const VirtualizedGridEntity& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    VirtualizedGridEntity& operator=( VirtualizedGridEntity&& ) = default;

    bool equals(const VirtualizedGridEntity& other) const
    {
      return impl_->equals(other);
    }

    //! returns true if father entity exists
    bool hasFather () const {
      return impl_->hasFather();
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return impl_->seed();
    }

    //! Level of this element
    int level () const
    {
      return impl_->level();
    }


    /** \brief The partition type for parallel computing */
    PartitionType partitionType () const {
      return impl_->partitionType();
    }


    //! Geometry of this entity
    Geometry geometry () const
    {
      return Geometry( impl_->geometry() );
    }


    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int codim) const
    {
      return impl_->subEntities(codim);
    }


    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template< int cc >
    typename GridImp::template Codim<cc>::Entity subEntity (int i) const {
      if constexpr (cc == 0)
        return impl_->subEntity0(i);
      if constexpr (cc == 1)
        return impl_->subEntity1(i);
      if constexpr (cc == dim-1)
        return impl_->subEntityDimMinus1(i);
      if constexpr (cc == dim)
        return impl_->subEntityDim(i);
    }

    //! First level intersection
    VirtualizedGridLevelIntersectionIterator<const GridImp> ilevelbegin () const {
      return VirtualizedGridLevelIntersectionIterator<const GridImp>(
        impl_->ilevelbegin());
    }


    //! Reference to one past the last neighbor
    VirtualizedGridLevelIntersectionIterator<const GridImp> ilevelend () const {
      return VirtualizedGridLevelIntersectionIterator<const GridImp>(
        impl_->ilevelend());
    }


    //! First leaf intersection
    VirtualizedGridLeafIntersectionIterator<const GridImp> ileafbegin () const {
      return VirtualizedGridLeafIntersectionIterator<const GridImp>(
        impl_->ileafbegin());
    }


    //! Reference to one past the last leaf intersection
    VirtualizedGridLeafIntersectionIterator<const GridImp> ileafend () const {
      return VirtualizedGridLeafIntersectionIterator<const GridImp>(
        impl_->ileafend());
    }


    //! returns true if Entity has NO children
    bool isLeaf() const {
      return impl_->isLeaf();
    }


    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    typename GridImp::template Codim<0>::Entity father () const {
      return impl_->father();
    }


    /** \brief Location of this element relative to the reference element element of the father.
     * This is sufficient to interpolate all dofs in conforming case.
     * Nonconforming may require access to neighbors of father and
     * computations with local coordinates.
     * On the fly case is somewhat inefficient since dofs  are visited several times.
     * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
     * implementation of numerical algorithms is only done for simple discretizations.
     * Assumes that meshes are nested.
     */
    LocalGeometry geometryInFather () const
    {
      return LocalGeometry( impl_->geometryInFather() );
    }


    /** \brief Inter-level access to son elements on higher levels<=maxlevel.
     * This is provided for sparsely stored nested unstructured meshes.
     * Returns iterator to first son.
     */
    VirtualizedGridHierarchicIterator<const GridImp> hbegin (int maxLevel) const
    {
      return VirtualizedGridHierarchicIterator<const GridImp>( impl_->hbegin(maxLevel) );
    }


    //! Returns iterator to one past the last son
    VirtualizedGridHierarchicIterator<const GridImp> hend (int maxLevel) const
    {
      return VirtualizedGridHierarchicIterator<const GridImp>( impl_->hend(maxLevel) );
    }

    std::unique_ptr<Interface> impl_;

  }; // end of VirtualizedGridEntity codim = 0


} // namespace Dune


#endif
