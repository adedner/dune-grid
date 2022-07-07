// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_VIRTUALIZEDGRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_HH

/** \file
 * \brief The VirtualizedGrid class
 */

#include <string>
#include <map>
#include <any>

#include <dune/common/parallel/communication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// The components of the VirtualizedGrid interface
#include "virtualizedgrid/cast.hh"
#include "virtualizedgrid/geometry.hh"
#include "virtualizedgrid/entity.hh"
#include "virtualizedgrid/entityseed.hh"
#include "virtualizedgrid/idtype.hh"
#include "virtualizedgrid/intersectioniterator.hh"
#include "virtualizedgrid/leveliterator.hh"
#include "virtualizedgrid/leafiterator.hh"
#include "virtualizedgrid/hierarchiciterator.hh"
#include "virtualizedgrid/indexsets.hh"

#if HAVE_MPI
  #include <dune/common/parallel/mpicommunication.hh>
  using VirtualizedCollectiveCommunication = Dune::Communication<MPI_Comm>;
#else
  using VirtualizedCollectiveCommunication = Dune::Communication<No_Comm>;
#endif

namespace Dune
{
  // Forward declaration
  template<int dimension, int dimensionworld, typename ct = double>
  class VirtualizedGrid;

  template<int dimension, int dimensionworld, typename ct>
  struct VirtualizedGridFamily
  {

  public:

    typedef GridTraits<
        dimension,
        dimensionworld,
        VirtualizedGrid<dimension, dimensionworld, ct>,
        VirtualizedGridGeometry,
        VirtualizedGridEntity,
        VirtualizedGridLevelIterator,
        VirtualizedGridLeafIntersection,
        VirtualizedGridLevelIntersection,
        VirtualizedGridLeafIntersectionIterator,
        VirtualizedGridLevelIntersectionIterator,
        VirtualizedGridHierarchicIterator,
        VirtualizedGridLeafIterator,
        VirtualizedGridLevelIndexSet<const VirtualizedGrid<dimension, dimensionworld, ct>>,
        VirtualizedGridLeafIndexSet<const VirtualizedGrid<dimension, dimensionworld, ct>>,
        VirtualizedGridGlobalIdSet<const VirtualizedGrid<dimension, dimensionworld, ct>>,
        VirtualizedGridIdType,
        VirtualizedGridLocalIdSet<const VirtualizedGrid<dimension, dimensionworld, ct>>,
        VirtualizedGridIdType,
        VirtualizedCollectiveCommunication,
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
  : public GridDefaultImplementation<dimension, dimensionworld, ct, VirtualizedGridFamily<dimension, dimensionworld, ct>>
  {
    typedef VirtualizedGrid<dimension, dimensionworld, ct> ThisType;

    friend class VirtualizedGridLevelIndexSet<const ThisType>;
    friend class VirtualizedGridLeafIndexSet<const ThisType>;
    friend class VirtualizedGridGlobalIdSet<const ThisType>;
    friend class VirtualizedGridLocalIdSet<const ThisType>;
    friend class VirtualizedGridHierarchicIterator<const ThisType>;
    friend class VirtualizedGridLevelIntersectionIterator<const ThisType>;
    friend class VirtualizedGridLeafIntersectionIterator<const ThisType>;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class VirtualizedGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class VirtualizedGridLeafIterator;

    template<int codim_, int dim_, class GridImp_>
    friend class VirtualizedGridEntity;

  public:
    //! type of the used GridFamily for this grid
    typedef VirtualizedGridFamily<dimension, dimensionworld, ct> GridFamily;

    //! the Traits
    typedef typename GridFamily::Traits Traits;

    //! The type used to store coordinates
    typedef ct ctype;

  private:
    typedef typename Traits::template Codim<0>::EntitySeed EntitySeed0;
    typedef typename Traits::template Codim<1>::EntitySeed EntitySeed1;
    typedef typename Traits::template Codim<dimension-1>::EntitySeed EntitySeedDimMinus1;
    typedef typename Traits::template Codim<dimension>::EntitySeed EntitySeedDim;

    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;

      virtual int maxLevel() const = 0;
      virtual typename Traits::template Codim<0>::LevelIterator lbegin (int level) const = 0;
      virtual typename Traits::template Codim<0>::LevelIterator lend (int level) const = 0;
      virtual typename Traits::template Codim<1>::LevelIterator lbegin1 (int level) const = 0;
      virtual typename Traits::template Codim<1>::LevelIterator lend1 (int level) const = 0;
      virtual typename Traits::template Codim<dimension-1>::LevelIterator lbeginDimMinus1 (int level) const = 0;
      virtual typename Traits::template Codim<dimension-1>::LevelIterator lendDimMinus1 (int level) const = 0;
      virtual typename Traits::template Codim<dimension>::LevelIterator lbeginDim (int level) const = 0;
      virtual typename Traits::template Codim<dimension>::LevelIterator lendDim (int level) const = 0;
      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LevelIterator lbeginGhost (int level) const = 0;
      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LevelIterator lendGhost (int level) const = 0;
      virtual typename Traits::template Codim<0>::template Partition<InteriorBorder_Partition>::LevelIterator lbeginInteriorBorder (int level) const = 0;
      virtual typename Traits::template Codim<0>::template Partition<InteriorBorder_Partition>::LevelIterator lendInteriorBorder (int level) const = 0;
      virtual typename Traits::template Codim<0>::LeafIterator leafbegin () const = 0;
      virtual typename Traits::template Codim<0>::LeafIterator leafend () const = 0;
      virtual typename Traits::template Codim<1>::LeafIterator leafbegin1 () const = 0;
      virtual typename Traits::template Codim<1>::LeafIterator leafend1 () const = 0;
      virtual typename Traits::template Codim<dimension-1>::LeafIterator leafbeginDimMinus1 () const = 0;
      virtual typename Traits::template Codim<dimension-1>::LeafIterator leafendDimMinus1 () const = 0;
      virtual typename Traits::template Codim<dimension>::LeafIterator leafbeginDim () const = 0;
      virtual typename Traits::template Codim<dimension>::LeafIterator leafendDim () const = 0;
      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LeafIterator leafbeginGhost() const = 0;
      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LeafIterator leafendGhost() const = 0;
      virtual typename Traits::template Codim<1>::template Partition<Ghost_Partition>::LeafIterator leafbegin1Ghost() const = 0;
      virtual typename Traits::template Codim<1>::template Partition<Ghost_Partition>::LeafIterator leafend1Ghost() const = 0;
      virtual int size (int level, int codim) const = 0;
      virtual size_t numBoundarySegments () const = 0;
      virtual int size (int codim) const = 0;
      virtual int size (int level, GeometryType type) const = 0;
      virtual int size (GeometryType type) const = 0;
      virtual const typename Traits::GlobalIdSet& globalIdSet() const = 0;
      virtual const typename Traits::LocalIdSet& localIdSet() const = 0;
      virtual const typename Traits::LevelIndexSet& levelIndexSet(int level) const = 0;
      virtual const typename Traits::LeafIndexSet& leafIndexSet() const = 0;
      virtual typename Traits::template Codim<0>::Entity entity0(const EntitySeed0& seed) const = 0;
      virtual typename Traits::template Codim<1>::Entity entity1(const EntitySeed1& seed) const = 0;
      virtual typename Traits::template Codim<dimension-1>::Entity entityDimMinus1(const EntitySeedDimMinus1& seed) const = 0;
      virtual typename Traits::template Codim<dimension>::Entity entityDim(const EntitySeedDim& seed) const = 0;
      virtual void globalRefine (int refCount) = 0;
      virtual bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) = 0;
      virtual int getMark(const typename Traits::template Codim<0>::Entity & e) const = 0;
      virtual bool preAdapt() = 0;
      virtual bool adapt() = 0;
      virtual void postAdapt() = 0;
      virtual unsigned int overlapSize(int codim) const = 0;
      virtual unsigned int ghostSize(int codim) const = 0;
      virtual unsigned int overlapSize(int level, int codim) const = 0;
      virtual unsigned int ghostSize(int level, int codim) const = 0;
      virtual const VirtualizedCollectiveCommunication& comm () const = 0;

    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename VirtualizedGridEntity<0, dimension, const ThisType>::template Implementation<const typename std::decay_t<I>::template Codim<0>::Entity> ImplEntity;
      typedef typename VirtualizedGridEntitySeed<0, const ThisType>::template Implementation<typename std::decay_t<I>::template Codim<0>::EntitySeed> ImplSeed0;
      typedef typename VirtualizedGridEntitySeed<1, const ThisType>::template Implementation<typename std::decay_t<I>::template Codim<1>::EntitySeed> ImplSeed1;
      typedef typename VirtualizedGridEntitySeed<dimension-1, const ThisType>::template Implementation<typename std::decay_t<I>::template Codim<dimension-1>::EntitySeed> ImplSeedDimMinus1;
      typedef typename VirtualizedGridEntitySeed<dimension, const ThisType>::template Implementation<typename std::decay_t<I>::template Codim<dimension>::EntitySeed> ImplSeedDim;

      Implementation ( I&& i )
      : impl_( std::forward<I>(i) ),
        globalIdSet_( impl().globalIdSet() ),
        localIdSet_( impl().localIdSet() ),
        leafIndexSet_( impl().leafIndexSet() )
      {
        for (int i = 0; i <= maxLevel(); i++)
        {
          VirtualizedGridLevelIndexSet<const ThisType>* p
            = new VirtualizedGridLevelIndexSet<const ThisType>( impl().levelIndexSet(i) );
          levelIndexSets_.push_back(p);
        }
      }

      ~Implementation()
      {
        for (size_t i = 0; i < levelIndexSets_.size(); i++)
          if (levelIndexSets_[i])
            delete (levelIndexSets_[i]);
      }

      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual int maxLevel() const override { return impl().maxLevel(); }

      virtual typename Traits::template Codim<0>::LevelIterator lbegin (int level) const override
      {
        return VirtualizedGridLevelIterator<0, All_Partition, const ThisType> ( std::move( impl().template lbegin<0>(level) ) );
      }

      virtual typename Traits::template Codim<0>::LevelIterator lend (int level) const override
      {
        return VirtualizedGridLevelIterator<0, All_Partition, const ThisType> ( std::move( impl().template lend<0>(level) ) );
      }

      virtual typename Traits::template Codim<1>::LevelIterator lbegin1 (int level) const override
      {
        return VirtualizedGridLevelIterator<1, All_Partition, const ThisType> ( std::move( impl().template lbegin<1>(level) ) );
      }

      virtual typename Traits::template Codim<1>::LevelIterator lend1 (int level) const override
      {
        return VirtualizedGridLevelIterator<1, All_Partition, const ThisType> ( std::move( impl().template lend<1>(level) ) );
      }

      virtual typename Traits::template Codim<dimension-1>::LevelIterator lbeginDimMinus1 (int level) const override
      {
        return VirtualizedGridLevelIterator<dimension-1, All_Partition, const ThisType> ( std::move( impl().template lbegin<dimension-1>(level) ) );
      }

      virtual typename Traits::template Codim<dimension-1>::LevelIterator lendDimMinus1 (int level) const override
      {
        return VirtualizedGridLevelIterator<dimension-1, All_Partition, const ThisType> ( std::move( impl().template lend<dimension-1>(level) ) );
      }

      virtual typename Traits::template Codim<dimension>::LevelIterator lbeginDim (int level) const override
      {
        return VirtualizedGridLevelIterator<dimension, All_Partition, const ThisType> ( std::move( impl().template lbegin<dimension>(level) ) );
      }

      virtual typename Traits::template Codim<dimension>::LevelIterator lendDim (int level) const override
      {
        return VirtualizedGridLevelIterator<dimension, All_Partition, const ThisType> ( std::move( impl().template lend<dimension>(level) ) );
      }

      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LevelIterator lbeginGhost (int level) const override
      {
        return VirtualizedGridLevelIterator<0, Ghost_Partition, const ThisType> ( std::move( impl().template lbegin<0, Ghost_Partition>(level) ) );
      }

      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LevelIterator lendGhost (int level) const override
      {
        return VirtualizedGridLevelIterator<0, Ghost_Partition, const ThisType> ( std::move( impl().template lend<0, Ghost_Partition>(level) ) );
      }

      virtual typename Traits::template Codim<0>::template Partition<InteriorBorder_Partition>::LevelIterator lbeginInteriorBorder (int level) const override
      {
        return VirtualizedGridLevelIterator<0, InteriorBorder_Partition, const ThisType> ( std::move( impl().template lbegin<0, InteriorBorder_Partition>(level) ) );
      }

      virtual typename Traits::template Codim<0>::template Partition<InteriorBorder_Partition>::LevelIterator lendInteriorBorder (int level) const override
      {
        return VirtualizedGridLevelIterator<0, InteriorBorder_Partition, const ThisType> ( std::move( impl().template lend<0, InteriorBorder_Partition>(level) ) );
      }

      virtual typename Traits::template Codim<0>::LeafIterator leafbegin () const override
      {
        return VirtualizedGridLeafIterator<0, All_Partition, const ThisType> ( std::move( impl().template leafbegin<0>() ) );
      }

      virtual typename Traits::template Codim<0>::LeafIterator leafend () const override
      {
        return VirtualizedGridLeafIterator<0, All_Partition, const ThisType> ( std::move( impl().template leafend<0>() ) );
      }

      virtual typename Traits::template Codim<1>::LeafIterator leafbegin1 () const override
      {
        return VirtualizedGridLeafIterator<1, All_Partition, const ThisType> ( std::move( impl().template leafbegin<1>() ) );
      }

      virtual typename Traits::template Codim<1>::LeafIterator leafend1 () const override
      {
        return VirtualizedGridLeafIterator<1, All_Partition, const ThisType> ( std::move( impl().template leafend<1>() ) );
      }

      virtual typename Traits::template Codim<dimension-1>::LeafIterator leafbeginDimMinus1 () const override
      {
        return VirtualizedGridLeafIterator<dimension-1, All_Partition, const ThisType> ( std::move( impl().template leafbegin<dimension-1>() ) );
      }

      virtual typename Traits::template Codim<dimension-1>::LeafIterator leafendDimMinus1 () const override
      {
        return VirtualizedGridLeafIterator<dimension-1, All_Partition, const ThisType> ( std::move( impl().template leafend<dimension-1>() ) );
      }

      virtual typename Traits::template Codim<dimension>::LeafIterator leafbeginDim () const override
      {
        return VirtualizedGridLeafIterator<dimension, All_Partition, const ThisType> ( std::move( impl().template leafbegin<dimension>() ) );
      }

      virtual typename Traits::template Codim<dimension>::LeafIterator leafendDim () const override
      {
        return VirtualizedGridLeafIterator<dimension, All_Partition, const ThisType> ( std::move( impl().template leafend<dimension>() ) );
      }

      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LeafIterator leafbeginGhost () const override
      {
        return VirtualizedGridLeafIterator<0, Ghost_Partition, const ThisType> ( std::move( impl().template leafbegin<0, Ghost_Partition>() ) );
      }

      virtual typename Traits::template Codim<0>::template Partition<Ghost_Partition>::LeafIterator leafendGhost () const override
      {
        return VirtualizedGridLeafIterator<0, Ghost_Partition, const ThisType> ( std::move( impl().template leafend<0, Ghost_Partition>() ) );
      }

      virtual typename Traits::template Codim<1>::template Partition<Ghost_Partition>::LeafIterator leafbegin1Ghost () const override
      {
        return VirtualizedGridLeafIterator<1, Ghost_Partition, const ThisType> ( std::move( impl().template leafbegin<1, Ghost_Partition>() ) );
      }

      virtual typename Traits::template Codim<1>::template Partition<Ghost_Partition>::LeafIterator leafend1Ghost () const override
      {
        return VirtualizedGridLeafIterator<1, Ghost_Partition, const ThisType> ( std::move( impl().template leafend<1, Ghost_Partition>() ) );
      }

      virtual int size (int level, int codim) const override { return impl().size(level, codim); }
      virtual size_t numBoundarySegments () const override { return impl().numBoundarySegments(); }
      virtual int size (int codim) const override { return impl().size(codim); }
      virtual int size (int level, GeometryType type) const override { return impl().size(level, type); }
      virtual int size (GeometryType type) const override { return impl().size(type); }

      virtual const typename Traits::GlobalIdSet& globalIdSet() const override
      {
        return dynamic_cast<const typename Traits::GlobalIdSet&>( globalIdSet_ );
      }

      virtual const typename Traits::LocalIdSet& localIdSet() const override
      {
        return dynamic_cast<const typename Traits::LocalIdSet&>( localIdSet_ );
      }

      virtual const typename Traits::LevelIndexSet& levelIndexSet(int level) const override
      {
        return dynamic_cast<const typename Traits::LevelIndexSet&>( *levelIndexSets_[level] );
      }

      virtual const typename Traits::LeafIndexSet& leafIndexSet() const override
      {
        return dynamic_cast<const typename Traits::LeafIndexSet&>( leafIndexSet_ );
      }

      virtual typename Traits::template Codim<0>::Entity entity0(const EntitySeed0& seed) const override
      {
        return VirtualizedGridEntity<0, dimension, const ThisType>( std::move( impl().entity(
          upcast<ImplSeed0>(seed)
        ) ) );
      }

      virtual typename Traits::template Codim<1>::Entity entity1(const EntitySeed1& seed) const override
      {
        return VirtualizedGridEntity<1, dimension, const ThisType>( std::move( impl().entity(
          upcast<ImplSeed1>(seed)
        ) ) );
      }

      virtual typename Traits::template Codim<dimension-1>::Entity entityDimMinus1(const EntitySeedDimMinus1& seed) const override
      {
        return VirtualizedGridEntity<dimension-1, dimension, const ThisType>( std::move( impl().entity(
          upcast<ImplSeedDimMinus1>(seed)
        ) ) );
      }

      virtual typename Traits::template Codim<dimension>::Entity entityDim(const EntitySeedDim& seed) const override
      {
        return VirtualizedGridEntity<dimension, dimension, const ThisType>( std::move( impl().entity(
          upcast<ImplSeedDim>(seed)
        ) ) );
      }

      virtual void globalRefine (int refCount) override { return impl().globalRefine(refCount); }
      virtual bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) override
      {
        return impl().mark(refCount,
          upcast<ImplEntity>(e)
        );
      }

      virtual int getMark(const typename Traits::template Codim<0>::Entity & e) const override
      {
        return impl().getMark(
          upcast<ImplEntity>(e)
        );
      }

      virtual bool preAdapt() override { return impl().preAdapt(); }
      virtual bool adapt() override { return impl().adapt(); }
      virtual void postAdapt() override { return impl().postAdapt(); }
      virtual unsigned int overlapSize(int codim) const override { return impl().overlapSize(codim); }
      virtual unsigned int ghostSize(int codim) const override { return impl().ghostSize(codim); }
      virtual unsigned int overlapSize(int level, int codim) const override { return impl().overlapSize(level, codim); }
      virtual unsigned int ghostSize(int level, int codim) const override { return impl().ghostSize(level, codim); }
      virtual const VirtualizedCollectiveCommunication& comm () const override { return impl().comm(); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I impl_;
      VirtualizedGridGlobalIdSet<const ThisType> globalIdSet_;
      VirtualizedGridLocalIdSet<const ThisType> localIdSet_;
      std::vector<VirtualizedGridLevelIndexSet<const ThisType>*> levelIndexSets_;
      VirtualizedGridLeafIndexSet<const ThisType> leafIndexSet_;
    };
    // VIRTUALIZATION END

  public:

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    /** \brief Constructor
     *
     * \param grid The grid hold by the VirtualizedGrid
     */
    template<class Impl>
    VirtualizedGrid(Impl&& grid)
    : impl_( new Implementation< Impl >( std::forward<Impl>(grid) ) )
    {}

    VirtualizedGrid(const VirtualizedGrid& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGrid ( VirtualizedGrid && ) = default;

    VirtualizedGrid& operator=(const VirtualizedGrid& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
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
      return impl_->lbegin(level);
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return impl_->lend(level);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      if constexpr (codim == 0)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lbegin(level);
        if constexpr (PiType == Ghost_Partition)
          return impl_->lbeginGhost(level);
        if constexpr (PiType == InteriorBorder_Partition)
          return impl_->lbeginInteriorBorder(level);
      }
      if constexpr (codim == 1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lbegin1(level);
      }
      if constexpr (codim == dimension-1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lbeginDimMinus1(level);
      }
      if constexpr (codim == dimension)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lbeginDim(level);
      }
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      if constexpr (codim == 0)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lend(level);
        if constexpr (PiType == Ghost_Partition)
          return impl_->lendGhost(level);
        if constexpr (PiType == InteriorBorder_Partition)
          return impl_->lendInteriorBorder(level);
      }
      if constexpr (codim == 1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lend1(level);
      }
      if constexpr (codim == dimension-1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lendDimMinus1(level);
      }
      if constexpr (codim == dimension)
      {
        if constexpr (PiType == All_Partition)
          return impl_->lendDim(level);
      }
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return impl_->leafbegin();
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return impl_->leafend();
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      if constexpr (codim == 0)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafbegin();
        if constexpr (PiType == Ghost_Partition)
          return impl_->leafbeginGhost();
      }
      if constexpr (codim == 1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafbegin1();
        if constexpr (PiType == Ghost_Partition)
          return impl_->leafbegin1Ghost();
      }
      if constexpr (codim == dimension-1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafbeginDimMinus1();
      }
      if constexpr (codim == dimension)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafbeginDim();
      }
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      if constexpr (codim == 0)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafend();
        if constexpr (PiType == Ghost_Partition)
          return impl_->leafendGhost();
      }
      if constexpr (codim == 1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafend1();
        if constexpr (PiType == Ghost_Partition)
          return impl_->leafend1Ghost();
      }
      if constexpr (codim == dimension-1)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafendDimMinus1();
      }
      if constexpr (codim == dimension)
      {
        if constexpr (PiType == All_Partition)
          return impl_->leafendDim();
      }
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
      return impl_->size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      return impl_->size(level, type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return impl_->size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const {
      return impl_->globalIdSet();
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const {
      return impl_->localIdSet();
    }


    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      return impl_->levelIndexSet(level);
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return impl_->leafIndexSet();
    }


    /** \brief Create Entity from EntitySeed */
    template < class EntitySeed >
    typename Traits::template Codim<EntitySeed::codimension>::Entity
    entity(const EntitySeed& seed) const
    {
      if constexpr (EntitySeed::codimension == 0)
        return impl_->entity0(seed);
      if constexpr (EntitySeed::codimension == 1)
        return impl_->entity1(seed);
      if constexpr (EntitySeed::codimension == dimension-1)
        return impl_->entityDimMinus1(seed);
      if constexpr (EntitySeed::codimension == dimension)
        return impl_->entityDim(seed);
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
      return this->leafGridView().overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return this->leafGridView().ghostSize(codim);
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return this->levelGridView(level).overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return this->levelGridView(level).ghostSize(codim);
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


    //! Returns the collective communication object
    const VirtualizedCollectiveCommunication& comm () const
    {
      return ccobj;
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype, CommunicationDirection dir) const
    {
      // TODO: we have to explicitly specify all CommDataHandleIF in virtual Interface class
      // impl_->communicate(data, iftype, dir);
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Returns the grid this VirtualizedGrid holds
    Interface& impl() const
    {
      return *impl_;
    }

  private:

    //! The grid this VirtualizedGrid holds
    std::unique_ptr< Interface > impl_;

    VirtualizedCollectiveCommunication ccobj;
  }; // end Class VirtualizedGrid




  namespace Capabilities
  {
    /** \brief has entities for all codimensions
     * \ingroup VirtualizedGrid
     */
    template<int dimension, int dimensionworld, typename ct, int codim>
    struct hasEntity<VirtualizedGrid<dimension, dimensionworld, ct>, codim>
    {
      static const bool v = true; // we assume this
    };

    template<int dimension, int dimensionworld, typename ct, int codim>
    struct hasEntityIterator<VirtualizedGrid<dimension, dimensionworld, ct>, codim>
    {
      static const bool v = true; // we assume this
    };

    /** \brief VirtualizedGrid can communicate
     *  \ingroup VirtualizedGrid
     */
    template<int dimension, int dimensionworld, typename ct, int codim>
    struct canCommunicate<VirtualizedGrid<dimension, dimensionworld, ct>, codim>
    {
      static const bool v = true; // we assume this
    };

    /** \brief has conforming level grids
     * \ingroup VirtualizedGrid
     */
    template<int dimension, int dimensionworld, typename ct>
    struct isLevelwiseConforming<VirtualizedGrid<dimension, dimensionworld, ct>>
    {
      static const bool v = true; // we assume this
    };

    /** \brief has conforming leaf grids
     * \ingroup VirtualizedGrid
     */
    template<int dimension, int dimensionworld, typename ct>
    struct isLeafwiseConforming<VirtualizedGrid<dimension, dimensionworld, ct>>
    {
      static const bool v = true; // we assume this
    };
  } // end namespace Capabilities

} // namespace Dune

#include "io/file/dgfparser/dgfvirtualized.hh"

#endif // DUNE_GRID_VIRTUALIZEDGRID_HH
