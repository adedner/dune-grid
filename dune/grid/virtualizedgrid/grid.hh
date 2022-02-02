// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_VIRTUALIZEDGRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_HH

/** \file
 * \brief The VirtualizedGrid class
 */

#include <string>
#include <map>

#include <dune/common/parallel/communication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// The components of the VirtualizedGrid interface
#include "virtualizedgrid/virtualizedgridgeometry.hh"
#include "virtualizedgrid/virtualizedgridentity.hh"
#include "virtualizedgrid/virtualizedgridentityseed.hh"
#include "virtualizedgrid/virtualizedgridintersectioniterator.hh"
#include "virtualizedgrid/virtualizedgridleveliterator.hh"
#include "virtualizedgrid/virtualizedgridleafiterator.hh"
#include "virtualizedgrid/virtualizedgridhierarchiciterator.hh"
#include "virtualizedgrid/virtualizedgridindexsets.hh"

namespace Dune
{
  // Forward declaration
  class VirtualizedGrid;

  template<int dimension, int dimensionworld>
  struct VirtualizedGridFamily
  {

  public:

    typedef GridTraits<
        dimension,
        dimensionworld,
        Dune::VirtualizedGrid,
        VirtualizedGridGeometry,
        VirtualizedGridEntity,
        VirtualizedGridLevelIterator,
        VirtualizedGridLeafIntersection,
        VirtualizedGridLevelIntersection,
        VirtualizedGridLeafIntersectionIterator,
        VirtualizedGridLevelIntersectionIterator,
        VirtualizedGridHierarchicIterator,
        VirtualizedGridLeafIterator,
        VirtualizedGridLevelIndexSet<const VirtualizedGrid>,
        VirtualizedGridLeafIndexSet<const VirtualizedGrid>,
        VirtualizedGridGlobalIdSet<const VirtualizedGrid>,
        std::size_t, // TODO: IdType
        VirtualizedGridLocalIdSet<const VirtualizedGrid>,
        std::size_t, // TODO: LocalIdType
        Communication<No_Comm>,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        VirtualizedGridEntitySeed
        > Traits;

  };

  //**********************************************************************
  //
  // --VirtualizedGrid
  //
  //************************************************************************
  /*!
   * \brief Provides a virtualized grid
   * \ingroup GridImplementations
   * \ingroup VirtualizedGrid
   */

  template<int dimension, int dimensionworld, typename ct>
  class VirtualizedGrid
  : public GridDefaultImplementation<dimension, dimensionworld, ct, VirtualizedGridFamily<dimension, dimensionworld>>
  {
    friend class VirtualizedGridLevelIndexSet<const VirtualizedGrid>;
    friend class VirtualizedGridLeafIndexSet<const VirtualizedGrid>;
    friend class VirtualizedGridGlobalIdSet<const VirtualizedGrid>;
    friend class VirtualizedGridLocalIdSet<const VirtualizedGrid>;
    friend class VirtualizedGridHierarchicIterator<const VirtualizedGrid>;
    friend class VirtualizedGridLevelIntersectionIterator<const VirtualizedGrid>;
    friend class VirtualizedGridLeafIntersectionIterator<const VirtualizedGrid>;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class VirtualizedGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class VirtualizedGridLeafIterator;


    template<int codim_, int dim_, class GridImp_>
    friend class VirtualizedGridEntity;

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;

      int maxLevel() const = 0;
      typename Traits::template Codim<0>::LevelIterator lbegin (int level) const = 0; // TODO: other codims
      typename Traits::template Codim<0>::LevelIterator lend (int level) const = 0; // TODO: other codims
      // typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const = 0; // TODO
      // typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const = 0; // TODO
      typename Traits::template Codim<0>::LeafIterator leafbegin() const = 0; // TODO: other codims
      typename Traits::template Codim<0>::LeafIterator leafend() const = 0; // TODO: other codims
      // typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const = 0; // TODO
      // typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const = 0; // TODO
      int size (int level, int codim) const = 0;
      size_t numBoundarySegments () const = 0;
      int size (int codim) const = 0;
      int size (int level, GeometryType type) const = 0;
      int size (GeometryType type) const = 0;
      const typename Traits::GlobalIdSet& globalIdSet() const = 0;
      const typename Traits::LocalIdSet& localIdSet() const = 0;
      const typename Traits::LevelIndexSet& levelIndexSet(int level) const = 0;
      const typename Traits::LeafIndexSet& leafIndexSet() const = 0;
      typename Traits::template Codim<0>::Entity entity(const EntitySeed& seed) const = 0; // TODO: other codims
      void globalRefine (int refCount) = 0;
      bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) = 0;
      int getMark(const typename Traits::template Codim<0>::Entity & e) const = 0;
      bool preAdapt() = 0;
      bool adapt() = 0;
      void postAdapt() = 0;
      unsigned int overlapSize(int codim) const = 0;
      unsigned int ghostSize(int codim) const = 0;
      unsigned int overlapSize(int level, int codim) const = 0;
      unsigned int ghostSize(int level, int codim) const = 0;
      const Communication<No_Comm>& comm () const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual Geometry geometry () const override { return impl().geometry(); }

      int maxLevel() const { return impl().geometry(); }
      typename Traits::template Codim<0>::LevelIterator lbegin (int level) const { return impl().lbegin(level); } // TODO: other codims
      typename Traits::template Codim<0>::LevelIterator lend (int level) const { return impl().lend(level); } // TODO: other codims
      // typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const { return impl().lbegin(level); } // TODO
      // typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const { return impl().lend(level); } // TODO
      typename Traits::template Codim<0>::LeafIterator leafbegin() const { return impl().leafbegin(); } // TODO: other codims
      typename Traits::template Codim<0>::LeafIterator leafend() const { return impl().leafend(); } // TODO: other codims
      // typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const { return impl().leafbegin(); } // TODO
      // typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const { return impl().leafend(); } // TODO
      int size (int level, int codim) const { return impl().size(level, codim); }
      size_t numBoundarySegments () const { return impl().numBoundarySegments(); }
      int size (int codim) const { return impl().size(codim); }
      int size (int level, GeometryType type) const { return impl().size(level, type); }
      int size (GeometryType type) const { return impl().size(type); }
      const typename Traits::GlobalIdSet& globalIdSet() const { return impl().globalIdSet(); }
      const typename Traits::LocalIdSet& localIdSet() const { return impl().localIdSet(); }
      const typename Traits::LevelIndexSet& levelIndexSet(int level) const { return impl().levelIndexSet(level); }
      const typename Traits::LeafIndexSet& leafIndexSet() const { return impl().leafIndexSet(); }
      typename Traits::template Codim<0>::Entity entity(const EntitySeed& seed) const { return impl().entity(seed); } // TODO: other codims
      void globalRefine (int refCount) { return impl().globalRefine(refCount); }
      bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) { return impl().mark(refCount, e); }
      int getMark(const typename Traits::template Codim<0>::Entity & e) const { return impl().getMark(e); }
      bool preAdapt() { return impl().preAdapt(); }
      bool adapt() { return impl().adapt(); }
      void postAdapt() { return impl().postAdapt(); }
      unsigned int overlapSize(int codim) const { return impl().overlapSize(codim); }
      unsigned int ghostSize(int codim) const { return impl().ghostSize(codim); }
      unsigned int overlapSize(int level, int codim) const { return impl().overlapSize(level, codim); }
      unsigned int ghostSize(int level, int codim) const { return impl().ghostSize(level, codim); }
      const Communication<No_Comm>& comm () const { return impl().comm(); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I& impl_;
    };
    // VIRTUALIZATION END

  public:

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! type of the used GridFamily for this grid
    typedef VirtualizedGridFamily<dimension, dimensionworld> GridFamily;

    //! the Traits
    typedef GridFamily::Traits Traits;

    //! The type used to store coordinates
    typedef ct ctype;

    /** \brief Constructor
     *
     * \param grid The grid hold by the VirtualizedGrid
     */
    template<class Impl>
    VirtualizedGrid(const Impl& grid)
    : impl_( new Implementation< Impl >( grid ) ),
      leafIndexSet_(*this),
      globalIdSet_(*this),
      localIdSet_(*this)
    {
      setIndices();
    }

    //! Desctructor
    ~VirtualizedGrid()
    {
      // Delete level index sets
      for (size_t i=0; i<levelIndexSets_.size(); i++)
        if (levelIndexSets_[i])
          delete (levelIndexSets_[i]);
    }


    /** \brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {
      return impl_->maxLevel();
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return VirtualizedGridLevelIterator<codim,All_Partition, const VirtualizedGrid>(this, level);
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return VirtualizedGridLevelIterator<codim,All_Partition, const VirtualizedGrid>(this, level, true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return VirtualizedGridLevelIterator<codim,PiType, const VirtualizedGrid>(this, level);
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return VirtualizedGridLevelIterator<codim,PiType, const VirtualizedGrid>(this, level, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return VirtualizedGridLeafIterator<codim,All_Partition, const VirtualizedGrid>(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return VirtualizedGridLeafIterator<codim,All_Partition, const VirtualizedGrid>(this, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return VirtualizedGridLeafIterator<codim,PiType, const VirtualizedGrid>(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return VirtualizedGridLeafIterator<codim,PiType, const VirtualizedGrid>(this, true);
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      return impl_->size(level,codim);
    }

    /** \brief returns the number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const {
      return impl_->numBoundarySegments();
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const {
      return leafIndexSet().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      return levelIndexSets_[level]->size(type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const {
      return globalIdSet_;
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const {
      return localIdSet_;
    }


    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level < 0 || level > maxLevel())
      {
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      }
      return *levelIndexSets_[level];
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }


    /** \brief Create Entity from EntitySeed */
    template < class EntitySeed >
    typename Traits::template Codim<EntitySeed::codimension>::Entity
    entity(const EntitySeed& seed) const
    {
      typedef VirtualizedGridEntity<
        EntitySeed::codimension,
        dimension,
        const typename Traits::Grid
        > EntityImp;

      return EntityImp(this, impl_->entity(seed));
    }


    /** @name Grid Refinement Methods */
    /*@{*/


    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount)
    {
      impl_->globalRefine(refCount);
    }

    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     * The parameter is currently ignored
     *
     * \return <ul>
     * <li> true, if marking was succesfull </li>
     * <li> false, if marking was not possible </li>
     * </ul>
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e)
    {
      return impl_->mark(refCount, e);
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const
    {
      return impl_->getMark(e);
    }

    /** \brief returns true, if at least one entity is marked for adaption */
    bool preAdapt() {
      return impl_->preAdapt();
    }


    //! Triggers the grid refinement process
    bool adapt()
    {
      return impl_->adapt();
    }

    /** \brief Clean up refinement markers */
    void postAdapt() {
      return impl_->postAdapt();
    }

    /*@}*/

    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return impl_->leafGridView().overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return impl_->leafGridView().ghostSize(codim);
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return impl_->levelGridView(level).overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return impl_->levelGridView(level).ghostSize(codim);
    }


#if 0
    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "VirtualizedGrid::loadBalance()");
    }
#endif


    /** \brief dummy collective communication */
    const Communication<No_Comm>& comm () const
    {
      return ccobj;
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Returns the grid this VirtualizedGrid holds
    Impl& impl() const
    {
      return *impl_;
    }

  private:

    //! compute the grid indices and ids
    void setIndices()
    {
      localIdSet_.update();

      globalIdSet_.update();

      // //////////////////////////////////////////
      //   Create the index sets
      // //////////////////////////////////////////
      for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
        VirtualizedGridLevelIndexSet<const VirtualizedGrid>* p
          = new VirtualizedGridLevelIndexSet<const VirtualizedGrid>();
        levelIndexSets_.push_back(p);
      }

      for (int i=0; i<=maxLevel(); i++)
        if (levelIndexSets_[i])
          levelIndexSets_[i]->update(*this, i);

      leafIndexSet_.update(*this);

    }

    //! The grid this VirtualizedGrid holds
    std::unique_ptr< Interface > impl_;

    Communication<No_Comm> ccobj;
    std::vector<VirtualizedGridLevelIndexSet<const VirtualizedGrid>*> levelIndexSets_;
    VirtualizedGridLeafIndexSet<const VirtualizedGrid> leafIndexSet_;
    VirtualizedGridGlobalIdSet<const VirtualizedGrid> globalIdSet_;
    VirtualizedGridLocalIdSet<const VirtualizedGrid> localIdSet_;

  }; // end Class VirtualizedGrid




  namespace Capabilities
  {
    /** \brief has entities for all codimensions
     * \ingroup VirtualizedGrid
     */
    template<int codim>
    struct hasEntity<VirtualizedGrid, codim>
    {
      static const bool v = true; // TODO
    };

    template<int codim>
    struct hasEntityIterator<VirtualizedGrid, codim>
    {
      static const bool v = true; // TODO
    };

    /** \brief VirtualizedGrid can communicate
     *  \ingroup VirtualizedGrid
     */
    template<int codim>
    struct canCommunicate<VirtualizedGrid, codim>
    {
      static const bool v = true; // TODO
    };

    /** \brief has conforming level grids
     * \ingroup VirtualizedGrid
     */
    struct isLevelwiseConforming<VirtualizedGrid>
    {
      static const bool v = true;
    };

    /** \brief has conforming leaf grids
     * \ingroup VirtualizedGrid
     */
    struct isLeafwiseConforming<VirtualizedGrid>
    {
      static const bool v = true;
    };
  } // end namespace Capabilities

} // namespace Dune

#endif // DUNE_GRID_VIRTUALIZEDGRID_HH
