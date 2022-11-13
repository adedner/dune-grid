// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_VIRTUALIZEDGRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_HH

/** \file
 * \brief The VirtualizedGrid class
 */

#include <string>
#include <map>
#include <any>

#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
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
#include "virtualizedgrid/gridview.hh"

#if HAVE_MPI
  #include <dune/common/parallel/mpicommunication.hh>
  using VirtualizedCommunication = Dune::Communication<MPI_Comm>;
#else
  using VirtualizedCommunication = Dune::Communication<No_Comm>;
#endif

namespace Dune
{
  template <PartitionIteratorType p>
  struct Partition
      : std::integral_constant<PartitionIteratorType,p>
  {};

  // Forward declaration
  template<int dimension, int dimensionworld, typename ct = double,
           typename GID = VirtualizedGridIdType, typename LID = GID>
  class VirtualizedGrid;

  template<int dimension, int dimensionworld, typename ct, typename GID, typename LID>
  struct VirtualizedGridFamily
  {
    struct Traits
    {
      /** \brief The type that implements the grid. */
      typedef VirtualizedGrid< dimension, dimensionworld, ct, GID, LID > Grid;

      /** \brief The type of the intersection at the leafs of the grid. */
      typedef Dune::Intersection< const Grid, VirtualizedGridLeafIntersection< const Grid > > LeafIntersection;
      /** \brief The type of the intersection at the levels of the grid. */
      typedef Dune::Intersection< const Grid, VirtualizedGridLevelIntersection< const Grid > > LevelIntersection;
      /** \brief The type of the intersection iterator at the leafs of the grid. */
      typedef Dune::IntersectionIterator< const Grid, VirtualizedGridLeafIntersectionIterator< const Grid >, VirtualizedGridLeafIntersection< const Grid > > LeafIntersectionIterator;
      /** \brief The type of the intersection iterator at the levels of the grid. */
      typedef Dune::IntersectionIterator< const Grid, VirtualizedGridLevelIntersectionIterator< const Grid >, VirtualizedGridLevelIntersection< const Grid > > LevelIntersectionIterator;

      /** \brief The type of the  hierarchic iterator. */
      typedef Dune::EntityIterator< 0, const Grid, VirtualizedGridHierarchicIterator< const Grid > > HierarchicIterator;

      /**
       * \brief Traits associated with a specific codim.
       * \tparam cd The codimension.
       */
      template <int cd>
      struct Codim
      {
      public:
        /** \brief The type of the geometry associated with the entity.*/
        typedef Dune::Geometry< dimension-cd, dimensionworld, const Grid, VirtualizedGridGeometry > Geometry;
        /** \brief The type of the local geometry associated with the entity.*/
        typedef Dune::Geometry< dimension-cd, dimension, const Grid, VirtualizedGridGeometry > LocalGeometry;
        /** \brief The type of the entity. */
        typedef Dune::Entity< cd, dimension, const Grid, VirtualizedGridEntity > Entity;

        /** \brief The type of the entity seed of this codim.*/
        typedef Dune::EntitySeed< const Grid, VirtualizedGridEntitySeed<cd, const Grid> > EntitySeed;

        /**
         * \brief Traits associated with a specific grid partition type.
         * \tparam pitype The type of the grid partition.
         */
        template <PartitionIteratorType pitype>
        struct Partition
        {
          /** \brief The type of the iterator over the level entities of this codim on this partition. */
          typedef Dune::EntityIterator< cd, const Grid, VirtualizedGridLevelIterator< cd, pitype, const Grid > > LevelIterator;
          /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
          typedef Dune::EntityIterator< cd, const Grid, VirtualizedGridLeafIterator< cd, pitype, const Grid > > LeafIterator;
        };

        /** \brief The type of the iterator over all leaf entities of this codim. */
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

        /** \brief The type of the entity pointer for entities of this codim.*/
        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

      private:
        friend class Dune::Entity< cd, dimension, const Grid, VirtualizedGridEntity >;
      };

      /** \brief The type of the leaf grid view. */
      typedef Dune::GridView< VirtualizedGridLeafViewTraits< const Grid > > LeafGridView;
      /** \brief The type of the level grid view. */
      typedef Dune::GridView< VirtualizedGridLevelViewTraits< const Grid > > LevelGridView;

      /** \brief The type of the level index set. */
      typedef IndexSet< const Grid, VirtualizedGridLevelIndexSet< const Grid > > LevelIndexSet;
      /** \brief The type of the leaf index set. */
      typedef IndexSet< const Grid, VirtualizedGridLeafIndexSet< const Grid > > LeafIndexSet;
      /** \brief The type of the global id set. */
      typedef IdSet< const Grid, VirtualizedGridGlobalIdSet< const Grid >, GID> GlobalIdSet;
      /** \brief The type of the local id set. */
      typedef IdSet< const Grid, VirtualizedGridLocalIdSet< const Grid >, LID> LocalIdSet;

      /** \brief The type of the collective communication. */
      typedef VirtualizedCommunication Communication;
    };
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

  template<int dimension, int dimensionworld, typename ct, typename GID, typename LID>
  class VirtualizedGrid
  : public GridDefaultImplementation<dimension, dimensionworld, ct, VirtualizedGridFamily<dimension, dimensionworld, ct,GID, LID>>
  {
    typedef VirtualizedGrid ThisType;

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
    typedef VirtualizedGridFamily<dimension, dimensionworld, ct, GID, LID> GridFamily;

    //! the Traits
    typedef typename GridFamily::Traits Traits;

    //! The type used to store coordinates
    typedef ct ctype;

  private:
    // VIRTUALIZATION BEGIN
    template<int codim, PartitionIteratorType pitype>
    struct InterfaceCodimPartition
    {
      virtual ~InterfaceCodimPartition () = default;

      using LevelIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator;
      using LeafIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator;

      virtual LevelIterator lbegin (Codim<codim>, Partition<pitype>, int level) const = 0;
      virtual LevelIterator lend (Codim<codim>, Partition<pitype>, int level) const = 0;
      virtual LeafIterator leafbegin (Codim<codim>, Partition<pitype>) const = 0;
      virtual LeafIterator leafend (Codim<codim>, Partition<pitype>) const = 0;
    };

    template<int codim, PartitionIteratorType... pitypes>
    struct InterfaceCodimImpl
        : virtual InterfaceCodimPartition<codim,pitypes>...
    {
      virtual ~InterfaceCodimImpl () = default;

      using LevelIterator = typename Traits::template Codim<codim>::LevelIterator;
      using LeafIterator = typename Traits::template Codim<codim>::LeafIterator;
      using Entity = typename Traits::template Codim<codim>::Entity;
      using EntitySeed = typename Traits::template Codim<codim>::EntitySeed;

      virtual LevelIterator lbegin (Codim<codim>, int level) const = 0;
      virtual LevelIterator lend (Codim<codim>, int level) const = 0;
      virtual LeafIterator leafbegin (Codim<codim>) const = 0;
      virtual LeafIterator leafend (Codim<codim>) const = 0;
      virtual Entity entity (Codim<codim>, const EntitySeed& seed) const = 0;

      using InterfaceCodimPartition<codim,pitypes>::lbegin...;
      using InterfaceCodimPartition<codim,pitypes>::lend...;
      using InterfaceCodimPartition<codim,pitypes>::leafbegin...;
      using InterfaceCodimPartition<codim,pitypes>::leafend...;
    };

    template<int codim>
    using InterfaceCodim = InterfaceCodimImpl<codim,
      Interior_Partition, InteriorBorder_Partition, Overlap_Partition, OverlapFront_Partition,
      All_Partition, Ghost_Partition>;

    template<int... codims>
    struct InterfaceImpl
        : virtual InterfaceCodim<codims>...
    {
      using Entity0 = typename Traits::template Codim<0>::Entity;

      virtual ~InterfaceImpl () = default;
      virtual InterfaceImpl *clone () const = 0;
      virtual int maxLevel() const = 0;
      virtual int size (int level, int codim) const = 0;
      virtual int size (int codim) const = 0;
      virtual int size (int level, GeometryType type) const = 0;
      virtual int size (GeometryType type) const = 0;
      virtual size_t numBoundarySegments () const = 0;
      virtual const typename Traits::GlobalIdSet& globalIdSet() const = 0;
      virtual const typename Traits::LocalIdSet& localIdSet() const = 0;
      virtual const typename Traits::LevelIndexSet& levelIndexSet(int level) const = 0;
      virtual const typename Traits::LeafIndexSet& leafIndexSet() const = 0;
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
      virtual const VirtualizedCommunication& comm () const = 0;

      virtual typename Traits::LevelIntersectionIterator ilevelbegin (const Entity0& entity) const = 0;
      virtual typename Traits::LevelIntersectionIterator ilevelend (const Entity0& entity) const = 0;
      virtual typename Traits::LeafIntersectionIterator ileafbegin (const Entity0& entity) const = 0;
      virtual typename Traits::LeafIntersectionIterator ileafend (const Entity0& entity) const = 0;

      using InterfaceCodim<codims>::lbegin...;
      using InterfaceCodim<codims>::lend...;
      using InterfaceCodim<codims>::leafbegin...;
      using InterfaceCodim<codims>::leafend...;
      using InterfaceCodim<codims>::entity...;
    };

    template<class Seq>
    struct Interface_t;
    template<int... codims>
    struct Interface_t<std::integer_sequence<int,codims...>> { using type = InterfaceImpl<codims...>; };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dimension+1>>::type;

    template<class Derived, class I, int codim, PartitionIteratorType pitype>
    struct DUNE_PRIVATE ImplementationCodimPartition
        : virtual InterfaceCodimPartition<codim,pitype>
    {
      using LevelIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator;
      using VirtLevelIterator = VirtualizedGridLevelIterator<codim,pitype,const ThisType>;
      using LeafIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator;
      using VirtLeafIterator = VirtualizedGridLeafIterator<codim,pitype,const ThisType>;

      virtual LevelIterator lbegin (Codim<codim>, Partition<pitype>, int level) const final {
        return VirtLevelIterator( std::move(derived().impl().levelGridView(level).template begin<codim,pitype>()) );
      }
      virtual LevelIterator lend (Codim<codim>, Partition<pitype>, int level) const final {
        return VirtLevelIterator( std::move(derived().impl().levelGridView(level).template end<codim,pitype>()) );
      }
      virtual LeafIterator leafbegin (Codim<codim>, Partition<pitype>) const final {
        return VirtLeafIterator( std::move(derived().impl().leafGridView().template begin<codim,pitype>()) );
      }
      virtual LeafIterator leafend (Codim<codim>, Partition<pitype>) const final {
        return VirtLeafIterator( std::move(derived().impl().leafGridView().template end<codim,pitype>()) );
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Derived, class I, int codim, PartitionIteratorType... pitypes>
    struct DUNE_PRIVATE ImplementationCodimImpl
        : virtual InterfaceCodimImpl<codim,pitypes...>
        , public ImplementationCodimPartition<Derived,I,codim,pitypes>...
    {
      using LevelIterator = typename Traits::template Codim<codim>::LevelIterator;
      using VirtLevelIterator = VirtualizedGridLevelIterator<codim,All_Partition,const ThisType>;
      using LeafIterator = typename Traits::template Codim<codim>::LeafIterator;
      using VirtLeafIterator = VirtualizedGridLeafIterator<codim,All_Partition,const ThisType>;
      using Entity = typename Traits::template Codim<codim>::Entity;
      using VirtEntity = VirtualizedGridEntity<codim,dimension,const ThisType>;
      using EntitySeed = typename Traits::template Codim<codim>::EntitySeed;
      using EntitySeedImpl = typename VirtualizedGridEntitySeed<codim,const ThisType>::template Implementation<typename std::decay_t<I>::template Codim<codim>::EntitySeed>;

      virtual LevelIterator lbegin (Codim<codim>, int level) const final {
        return VirtLevelIterator( std::move(derived().impl().levelGridView(level).template begin<codim>()) );
      }
      virtual LevelIterator lend (Codim<codim>, int level) const final {
        return VirtLevelIterator( std::move(derived().impl().levelGridView(level).template end<codim>()) );
      }
      virtual LeafIterator leafbegin (Codim<codim>) const final {
        return VirtLeafIterator( std::move(derived().impl().leafGridView().template begin<codim>()) );
      }
      virtual LeafIterator leafend (Codim<codim>) const final {
        return VirtLeafIterator( std::move(derived().impl().leafGridView().template end<codim>()) );
      }
      virtual Entity entity (Codim<codim>, const EntitySeed& seed) const final {
        return VirtEntity( std::move(derived().impl().entity(upcast<EntitySeedImpl>(seed))) );
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Derived, class I, int codim>
    using ImplementationCodim = ImplementationCodimImpl<Derived,I,codim,
      Interior_Partition, InteriorBorder_Partition, Overlap_Partition, OverlapFront_Partition,
      All_Partition, Ghost_Partition>;

    template<class I, int... codims>
    struct DUNE_PRIVATE ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      using Entity0 = typename Traits::template Codim<0>::Entity;
      using Entity0Impl = typename VirtualizedGridEntity<0,dimension,const ThisType>::template Implementation<const typename std::decay_t<I>::template Codim<0>::Entity>;
      using VirtLevelIntersectionIterator = VirtualizedGridLevelIntersectionIterator<const ThisType>;
      using VirtLeafIntersectionIterator = VirtualizedGridLeafIntersectionIterator<const ThisType>;

      ImplementationImpl ( I&& i )
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

      ~ImplementationImpl ()
      {
        for (size_t i = 0; i < levelIndexSets_.size(); i++)
          if (levelIndexSets_[i])
            delete (levelIndexSets_[i]);
      }

      virtual ImplementationImpl *clone () const override { return new ImplementationImpl( *this ); }
      virtual int maxLevel () const override { return impl().maxLevel(); }

      virtual typename Traits::LevelIntersectionIterator ilevelbegin (const Entity0& entity) const override {
        return VirtLevelIntersectionIterator( std::move(impl().levelGridView(entity.level()).ibegin(upcast<Entity0Impl>(entity))) );
      }
      virtual typename Traits::LevelIntersectionIterator ilevelend (const Entity0& entity) const override {
        return VirtLevelIntersectionIterator( std::move(impl().levelGridView(entity.level()).iend(upcast<Entity0Impl>(entity))) );
      }
      virtual typename Traits::LeafIntersectionIterator ileafbegin (const Entity0& entity) const override {
        return VirtLeafIntersectionIterator( std::move(impl().leafGridView().ibegin(upcast<Entity0Impl>(entity))) );
      }
      virtual typename Traits::LeafIntersectionIterator ileafend (const Entity0& entity) const override {
        return VirtLeafIntersectionIterator( std::move(impl().leafGridView().iend(upcast<Entity0Impl>(entity))) );
      }

      virtual int size (int level, int codim) const override { return impl().size(level, codim); }
      virtual int size (int codim) const override { return impl().size(codim); }
      virtual int size (int level, GeometryType type) const override { return impl().size(level, type); }
      virtual int size (GeometryType type) const override { return impl().size(type); }
      virtual size_t numBoundarySegments () const override { return impl().numBoundarySegments(); }

      virtual const typename Traits::GlobalIdSet& globalIdSet () const override {
        return dynamic_cast<const typename Traits::GlobalIdSet&>(globalIdSet_);
      }
      virtual const typename Traits::LocalIdSet& localIdSet () const override {
        return dynamic_cast<const typename Traits::LocalIdSet&>(localIdSet_);
      }
      virtual const typename Traits::LevelIndexSet& levelIndexSet (int level) const override {
        return dynamic_cast<const typename Traits::LevelIndexSet&>(*levelIndexSets_[level]);
      }
      virtual const typename Traits::LeafIndexSet& leafIndexSet () const override {
        return dynamic_cast<const typename Traits::LeafIndexSet&>(leafIndexSet_);
      }

      virtual void globalRefine (int refCount) override { return impl().globalRefine(refCount); }
      virtual bool mark (int refCount, const Entity0& e) override {
        return impl().mark(refCount, upcast<Entity0Impl>(e));
      }
      virtual int getMark (const Entity0 & e) const override {
        return impl().getMark(upcast<Entity0Impl>(e));
      }

      virtual bool preAdapt () override { return impl().preAdapt(); }
      virtual bool adapt () override { return impl().adapt(); }
      virtual void postAdapt () override { return impl().postAdapt(); }
      virtual unsigned int overlapSize (int codim) const override { return impl().overlapSize(codim); }
      virtual unsigned int ghostSize (int codim) const override { return impl().ghostSize(codim); }
      virtual unsigned int overlapSize (int level, int codim) const override { return impl().overlapSize(level, codim); }
      virtual unsigned int ghostSize (int level, int codim) const override { return impl().ghostSize(level, codim); }
      virtual const VirtualizedCommunication& comm () const override { return impl().comm(); }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
      VirtualizedGridGlobalIdSet<const ThisType> globalIdSet_;
      VirtualizedGridLocalIdSet<const ThisType> localIdSet_;
      std::vector<VirtualizedGridLevelIndexSet<const ThisType>*> levelIndexSets_;
      VirtualizedGridLeafIndexSet<const ThisType> leafIndexSet_;
    };

    template<class I, class Seq>
    struct Implementation_t;
    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>> { using type = ImplementationImpl<I,codims...>; };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dimension+1>>::type;

    // VIRTUALIZATION END

    struct Cache
    {
      template<class Grid>
      void update (const Grid& grid)
      {
        maxLevel_ = grid.maxLevel();

        levelSizes_.resize(maxLevel_+1);
        for (int level = 0; level <= maxLevel_; ++level)
          for (int codim = 0; codim <= dimension; ++codim)
            levelSizes_[level][codim] = grid.size(level, codim);

        for (int codim = 0; codim <= dimension; ++codim)
          sizes_[codim] = grid.size(codim);

        numBoundarySegments_ = grid.numBoundarySegments();
      }

      //! maximum level defined in this grid
      int maxLevel () const {
        return maxLevel_;
      }

      //! number of grid entities per level and codim
      int size (int level, int codim) const {
        return levelSizes_[level][codim];
      }

      //! number of leaf entities per codim in this process
      int size (int codim) const {
        return sizes_[codim];
      }

      //! number of boundary segments within the macro grid
      std::size_t numBoundarySegments () const {
        return numBoundarySegments_;
      }

    private:
      int maxLevel_{};
      std::vector<std::array<std::size_t,dimension+1>> levelSizes_{};
      std::array<std::size_t,dimension+1> sizes_{};
      std::size_t numBoundarySegments_{};
    };


  public:

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    /** \brief Constructor
     *
     * \param grid The grid hold by the VirtualizedGrid
     */
    template<class Impl>
    VirtualizedGrid (Impl&& grid)
    : impl_( new Implementation< Impl >( std::forward<Impl>(grid) ) )
    {
      update(grid);
    }

    VirtualizedGrid (const VirtualizedGrid& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    , cache_( other.cache_ )
    {}

    VirtualizedGrid ( VirtualizedGrid && ) = default;

    VirtualizedGrid& operator= (const VirtualizedGrid& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      cache_ = other.cache_;
      return *this;
    }

    /** \brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel () const {
      return cache().maxLevel();
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return impl_->lbegin(Codim<codim>{},level);
    }

    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return impl_->lend(Codim<codim>{},level);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lbegin (int level) const {
      return impl_->lbegin(Codim<codim>{},Partition<pitype>{},level);
    }

    //! one past the end on this level
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lend (int level) const {
      return impl_->lend(Codim<codim>{},Partition<pitype>{},level);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin () const {
      return impl_->leafbegin(Codim<codim>{});
    }

    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend () const {
      return impl_->leafend(Codim<codim>{});
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafbegin () const {
      return impl_->leafbegin(Codim<codim>{},Partition<pitype>{});
    }

    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafend () const {
      return impl_->leafend(Codim<codim>{},Partition<pitype>{});
    }


    virtual typename Traits::LevelIntersectionIterator ilevelbegin (const typename Traits::template Codim<0>::Entity& entity) const {
      return impl_->ilevelbegin( entity );
    }

    virtual typename Traits::LevelIntersectionIterator ilevelend (const typename Traits::template Codim<0>::Entity& entity) const {
      return impl_->ilevelend( entity );
    }

    virtual typename Traits::LeafIntersectionIterator ileafbegin (const typename Traits::template Codim<0>::Entity& entity) const {
      return impl_->ileafbegin( entity );
    }

    virtual typename Traits::LeafIntersectionIterator ileafend (const typename Traits::template Codim<0>::Entity& entity) const {
      return impl_->ileafend( entity );
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      return cache().size(level, codim);
    }

    /** \brief returns the number of boundary segments within the macro grid
     */
    std::size_t numBoundarySegments () const {
      return cache().numBoundarySegments();
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const {
      return cache().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      return impl_->size(level, type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const {
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
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const {
      return impl_->levelIndexSet(level);
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const {
      return impl_->leafIndexSet();
    }


    /** \brief Create Entity from EntitySeed */
    template<class EntitySeed>
    typename Traits::template Codim<EntitySeed::codimension>::Entity entity(const EntitySeed& seed) const {
      return impl_->entity(Codim<EntitySeed::codimension>{}, seed);
    }


    /** @name Grid Refinement Methods */
    /*@{*/


    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount) {
      impl_->globalRefine(refCount);
      update(*impl_);
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
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e) {
      return impl_->mark(refCount, e);
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const {
      return impl_->getMark(e);
    }

    /** \brief returns true, if at least one entity is marked for adaption */
    bool preAdapt() {
      return impl_->preAdapt();
    }


    //! Triggers the grid refinement process
    bool adapt() {
      bool result = impl_->adapt();
      update(*impl_);
      return result;
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
    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement) {
      DUNE_THROW(NotImplemented, "VirtualizedGrid::loadBalance()");
      update(*impl_);
    }
#endif


    //! Returns the collective communication object
    const VirtualizedCommunication& comm () const {
      return ccobj;
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype, CommunicationDirection dir, int level) const {
      // TODO: we have to explicitly specify all CommDataHandleIF in virtual Interface class
      // impl_->communicate(data, iftype, dir);
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype, CommunicationDirection dir) const {
      // TODO: we have to explicitly specify all CommDataHandleIF in virtual Interface class
      // impl_->communicate(data, iftype, dir);
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Returns the grid this VirtualizedGrid holds
    Interface& impl() const {
      return *impl_;
    }

    const auto& cache() const {
#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
      return cache_;
#else
      return *impl_;
#endif
    }

    template<class Grid>
    void update (const Grid& grid) {
#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
      cache_.update(grid);
#endif
    }

  private:
    //! The grid this VirtualizedGrid holds
    std::unique_ptr< Interface > impl_;

    // some cached information
    Cache cache_;

    VirtualizedCommunication ccobj;
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
