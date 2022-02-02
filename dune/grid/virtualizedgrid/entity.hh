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
    public EntityDefaultImplementation <codim,dim,GridImp,VirtualizedGridEntity>
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
      virtual bool equals(const VirtualizedGridEntity& other) const = 0;
      virtual bool hasFather () const = 0;
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
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool equals(const VirtualizedGridEntity& other) const override { return impl().equals(other); }
      virtual bool hasFather () const override { return impl().hasFather(); }
      virtual EntitySeed seed () const override { return impl().seed(); }
      virtual int level () const override { return impl().level(); }
      virtual PartitionType partitionType () const override { return impl().partitionType(); }
      virtual unsigned int subEntities (unsigned int cc) const override { return impl().subEntities(cc); }
      virtual Geometry geometry () const override { return impl().geometry(); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END

  public:

    VirtualizedGridEntity()
      : virtualizedGrid_(nullptr)
    {}

    template< class ImplGridEntity >
    VirtualizedGridEntity(const GridImp* virtualizedGrid, const ImplGridEntity& implEntity)
      : virtualizedGrid_(virtualizedGrid),
        impl_( new Implementation(implEntity) )
    {}

    VirtualizedGridEntity(const VirtualizedGridEntity& original)
      : virtualizedGrid_(original.virtualizedGrid_),
        impl_(original.impl_)
    {}

    VirtualizedGridEntity& operator=(const VirtualizedGridEntity& original)
    {
      if (this != &original)
      {
        virtualizedGrid_ = original.virtualizedGrid_;
        impl_ = original.impl_;
      }
      return *this;
    }

    bool equals(const VirtualizedGridEntity& other) const
    {
      return impl_ == other.impl_;
    }

    //! returns true if father entity exists
    bool hasFather () const {
      return impl_->hasFather();
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return EntitySeed(*impl_);
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

  private:
    const GridImp* virtualizedGrid_;

  public:
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
  class VirtualizedGridEntity<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, VirtualizedGridEntity>
  {
  public:
    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on this level
    typedef VirtualizedGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;

    //! The Iterator over intersections on the leaf level
    typedef VirtualizedGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef VirtualizedGridHierarchicIterator<GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

  private:

    typedef typename GridImp::ctype ctype;

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool equals(const VirtualizedGridEntity& other) const = 0;
      virtual bool hasFather () const = 0;
      virtual EntitySeed seed () const = 0;
      virtual int level () const = 0;
      virtual PartitionType partitionType () const = 0;
      virtual Geometry geometry () const = 0;
      virtual unsigned int subEntities (unsigned int cc) const = 0;
      virtual typename GridImp::template Codim<1>::Entity subEntity (int i) const = 0; // TODO: other codims
      virtual LevelIntersectionIterator ilevelbegin () const = 0;
      virtual LevelIntersectionIterator ilevelend () const = 0;
      virtual LeafIntersectionIterator ileafbegin () const = 0;
      virtual LeafIntersectionIterator ileafend () const = 0;
      virtual bool isLeaf() const = 0;
      virtual typename GridImp::template Codim<0>::Entity father () const = 0;
      virtual LocalGeometry geometryInFather () const = 0;
      virtual HierarchicIterator hbegin (int maxLevel) const = 0;
      virtual HierarchicIterator hend (int maxLevel) const = 0;
      virtual bool wasRefined () const = 0;
      virtual bool mightBeCoarsened () const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual bool equals(const VirtualizedGridEntity& other) const override { return impl().equals(other); }
      virtual bool hasFather () const override { return impl().hasFather(); }
      virtual EntitySeed seed () const override { return impl().seed(); }
      virtual int level () const override { return impl().level(); }
      virtual PartitionType partitionType () const override { return impl().partitionType(); }
      virtual Geometry geometry () const override { return impl().geometry(); }
      virtual unsigned int subEntities (unsigned int cc) const override { return impl().subEntities(cc); }
      virtual typename GridImp::template Codim<1>::Entity subEntity (int i) const { return impl().subEntity(i); } // TODO: other codims
      virtual LevelIntersectionIterator ilevelbegin () const { return impl().ilevelbegin(); };
      virtual LevelIntersectionIterator ilevelend () const { return impl().ilevelend(); };
      virtual LeafIntersectionIterator ileafbegin () const { return impl().ileafbegin(); }
      virtual LeafIntersectionIterator ileafend () const { return impl().ileafend(); }
      virtual bool isLeaf() const { return impl().isLeaf(); }
      virtual typename GridImp::template Codim<0>::Entity father () const { return impl().father(); }
      virtual LocalGeometry geometryInFather () const { return impl().geometryInFather(); }
      virtual HierarchicIterator hbegin (int maxLevel) const { return impl().hbegin(maxLevel); }
      virtual HierarchicIterator hend (int maxLevel) const { return impl().hend(maxLevel); }
      virtual bool wasRefined () const { return impl().wasRefined(); }
      virtual bool mightBeCoarsened () const { return impl().mightBeCoarsened(); }
    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END

  public:

    VirtualizedGridEntity()
      : virtualizedGrid_(nullptr)
    {}

    template< class ImplGridEntity >
    VirtualizedGridEntity(const GridImp* virtualizedGrid, const ImplGridEntity& implEntity)
    : virtualizedGrid_(virtualizedGrid),
      impl_( new Implementation(implEntity) )
    {}

    VirtualizedGridEntity(const VirtualizedGridEntity& original)
    : virtualizedGrid_(original.virtualizedGrid_),
      impl_(original.impl_)
    {}

    VirtualizedGridEntity& operator=(const VirtualizedGridEntity& original)
    {
      if (this != &original)
      {
        virtualizedGrid_ = original.virtualizedGrid_;
        impl_ = original.impl_;
      }
      return *this;
    }

    bool equals(const VirtualizedGridEntity& other) const
    {
      return impl_ == other.impl_;
    }

    //! returns true if father entity exists
    bool hasFather () const {
      return impl_->hasFather();
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return EntitySeed(impl_);
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
    template<int cc>
    typename GridImp::template Codim<cc>::Entity subEntity (int i) const {
      return VirtualizedGridEntity<cc,dim,GridImp>(virtualizedGrid_, impl_->template subEntity<cc>(i));
    }


    //! First level intersection
    VirtualizedGridLevelIntersectionIterator<GridImp> ilevelbegin () const {
      return VirtualizedGridLevelIntersectionIterator<GridImp>(
        virtualizedGrid_,
        virtualizedGrid_->impl().levelGridView(level()).ibegin(impl_));
    }


    //! Reference to one past the last neighbor
    VirtualizedGridLevelIntersectionIterator<GridImp> ilevelend () const {
      return VirtualizedGridLevelIntersectionIterator<GridImp>(
        virtualizedGrid_,
        virtualizedGrid_->impl().levelGridView(level()).iend(impl_));
    }


    //! First leaf intersection
    VirtualizedGridLeafIntersectionIterator<GridImp> ileafbegin () const {
      return VirtualizedGridLeafIntersectionIterator<GridImp>(
        virtualizedGrid_,
        virtualizedGrid_->impl().leafGridView().ibegin(impl_));
    }


    //! Reference to one past the last leaf intersection
    VirtualizedGridLeafIntersectionIterator<GridImp> ileafend () const {
      return VirtualizedGridLeafIntersectionIterator<GridImp>(
        virtualizedGrid_,
        virtualizedGrid_->impl().leafGridView().iend(impl_));
    }


    //! returns true if Entity has NO children
    bool isLeaf() const {
      return impl_->isLeaf();
    }


    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    typename GridImp::template Codim<0>::Entity father () const {
      return VirtualizedGridEntity(virtualizedGrid_, impl_->father());
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
    VirtualizedGridHierarchicIterator<GridImp> hbegin (int maxLevel) const
    {
      return VirtualizedGridHierarchicIterator<const GridImp>(virtualizedGrid_, *this, maxLevel);
    }


    //! Returns iterator to one past the last son
    VirtualizedGridHierarchicIterator<GridImp> hend (int maxLevel) const
    {
      return VirtualizedGridHierarchicIterator<const GridImp>(virtualizedGrid_, *this, maxLevel, true);
    }


    //! \todo Please doc me !
    bool wasRefined () const
    {
      if (virtualizedGrid_->adaptationStep!=GridImp::adaptDone)
        return false;

      int level = this->level();
      int index = virtualizedGrid_->levelIndexSet(level).index(*this);
      return virtualizedGrid_->refinementMark_[level][index];
    }


    //! \todo Please doc me !
    bool mightBeCoarsened () const
    {
      return true;
    }


    // /////////////////////////////////////////
    //   Internal stuff
    // /////////////////////////////////////////

    const GridImp* virtualizedGrid_;
    std::unique_ptr<Interface> impl_;

  }; // end of VirtualizedGridEntity codim = 0


} // namespace Dune


#endif
