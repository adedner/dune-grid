// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_VIRTUALIZEDGRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_HH

/** \file
 * \brief The Polymorphic class
 */

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// The components of the Polymorphic interface
#include "polymorphicgrid/grid.hh"

namespace Dune
{
  // Forward declaration
  template<int dimension, int dimensionworld, typename ct>
  class PolymorphicGrid;

  template<int dimension, int dimensionworld, typename ct>
  struct PolymorphicGridFamily
  {
    struct Traits
    {
      /** \brief The type that implements the grid. */
      using Grid = PolymorphicGrid<dimension, dimensionworld, ct>;

      using IntersectionImp = PolymorphicIntersection<const Grid>;
      using IntersectionIteratorImp = PolymorphicIntersectionIterator<const Grid>;

      /** \brief The type of the intersection at the leafs of the grid. */
      using LeafIntersection = Dune::Intersection<const Grid, IntersectionImp>;

      /** \brief The type of the intersection at the levels of the grid. */
      using LevelIntersection = Dune::Intersection<const Grid, IntersectionImp>;

      /** \brief The type of the intersection iterator at the leafs of the grid. */
      using LeafIntersectionIterator = Dune::IntersectionIterator<const Grid, IntersectionIteratorImp, IntersectionImp>;

      /** \brief The type of the intersection iterator at the levels of the grid. */
      using LevelIntersectionIterator = Dune::IntersectionIterator<const Grid, IntersectionIteratorImp, IntersectionImp>;

      /** \brief The type of the  hierarchic iterator. */
      using HierarchicIterator = Dune::EntityIterator<0, const Grid, PolymorphicEntityIterator<0,const Grid> >;

      /**
       * \brief Traits associated with a specific codim.
       * \tparam cd The codimension.
       */
      template <int cd>
      struct Codim
      {
      public:
        /** \brief The type of the geometry associated with the entity.*/
        using Geometry = Dune::Geometry<dimension-cd, dimensionworld, const Grid, PolymorphicGeometry>;

        /** \brief The type of the local geometry associated with the entity.*/
        using LocalGeometry = Dune::Geometry<dimension-cd, dimension, const Grid, PolymorphicGeometry>;

        /** \brief The type of the entity. */
        using Entity = Dune::Entity<cd, dimension, const Grid, PolymorphicEntity>;

        /** \brief The type of the entity seed of this codim.*/
        using EntitySeed = Dune::EntitySeed<const Grid, PolymorphicEntitySeed<cd, const Grid> >;

        /**
         * \brief Traits associated with a specific grid partition type.
         * \tparam pitype The type of the grid partition.
         */
        template <PartitionIteratorType pitype>
        struct Partition
        {
          /** \brief The type of the iterator over the level entities of this codim on this partition. */
          using LevelIterator = Dune::EntityIterator<cd, const Grid, PolymorphicEntityIterator<cd, const Grid> >;

          /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
          using LeafIterator = Dune::EntityIterator<cd, const Grid, PolymorphicEntityIterator<cd, const Grid> >;
        };

        /** \brief The type of the iterator over all leaf entities of this codim. */
        using LeafIterator = typename Partition< All_Partition >::LeafIterator;

        /** \brief The type of the entity pointer for entities of this codim.*/
        using LevelIterator = typename Partition< All_Partition >::LevelIterator;

      private:
        friend class Dune::Entity<cd, dimension, const Grid, PolymorphicEntity>;
      };

      /** \brief The type of the leaf grid view. */
      using LeafGridView = Dune::GridView<PolymorphicLeafViewTraits<const Grid> >;

      /** \brief The type of the level grid view. */
      using LevelGridView = Dune::GridView<PolymorphicLevelViewTraits<const Grid> >;

      /** \brief The type of the level index set. */
      using LevelIndexSet = IndexSet<const Grid, PolymorphicIndexSet<const Grid> >;

      /** \brief The type of the leaf index set. */
      using LeafIndexSet = IndexSet<const Grid, PolymorphicIndexSet<const Grid> >;

      /** \brief The type of the global id set. */
      using GlobalIdSet = IdSet<const Grid, PolymorphicIdSet<const Grid>, PolymorphicIdType>;

      /** \brief The type of the local id set. */
      using LocalIdSet = IdSet<const Grid, PolymorphicIdSet<const Grid>, PolymorphicIdType>;

      /** \brief The type of the collective communication. */
      using Communication = PolymorphicCommunication;
    };
  };


  template<int dimension, int dimensionworld, typename ct>
  class PolymorphicGrid
    : public GridDefaultImplementation<dimension, dimensionworld, ct, PolymorphicGridFamily<dimension, dimensionworld, ct>>
  {
    using Self = PolymorphicGrid<dimension, dimensionworld, ct>;
    using Definition = PolymorphicGridDefinition<Self,dimension,dimensionworld,ct>;

    template <PartitionIteratorType p>
    using _Partition = typename Definition::template _Partition<p>;

    template <class I>
    using Implementation = typename Definition::template Implementation<I>;

    using Interface = typename Definition::Interface;

  public:
    //! type of the used GridFamily for this grid
    using GridFamily = PolymorphicGridFamily<dimension, dimensionworld, ct>;

    //! the Traits
    using Traits = typename GridFamily::Traits;

    //! The type used to store coordinates
    using ctype = ct;

  public:

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    /** \brief Constructor
     *
     * \param grid The grid hold by the Polymorphic
     */
    template<class Impl>
    PolymorphicGrid (Impl&& grid)
      : impl_(new Implementation<Impl>( std::forward<Impl>(grid) ) )
    {}

    PolymorphicGrid (const PolymorphicGrid& other)
      : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    PolymorphicGrid (PolymorphicGrid &&) = default;

    PolymorphicGrid& operator= (const PolymorphicGrid& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }


    /** \brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const
    {
      return impl().maxLevel();
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const
    {
      return impl().lbegin(Codim<codim>{},level);
    }

    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const
    {
      return impl().lend(Codim<codim>{},level);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      return impl().lbegin(Codim<codim>{},_Partition<pitype>{},level);
    }

    //! one past the end on this level
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      return impl().lend(Codim<codim>{},_Partition<pitype>{},level);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin () const
    {
      return impl().leafbegin(Codim<codim>{});
    }

    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend () const
    {
      return impl().leafend(Codim<codim>{});
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      return impl().leafbegin(Codim<codim>{},_Partition<pitype>{});
    }

    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafend () const
    {
      return impl().leafend(Codim<codim>{},_Partition<pitype>{});
    }


    virtual typename Traits::LevelIntersectionIterator ilevelbegin (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return impl().ilevelbegin( entity );
    }

    virtual typename Traits::LevelIntersectionIterator ilevelend (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return impl().ilevelend( entity );
    }

    virtual typename Traits::LeafIntersectionIterator ileafbegin (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return impl().ileafbegin( entity );
    }

    virtual typename Traits::LeafIntersectionIterator ileafend (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return impl().ileafend( entity );
    }


    /** \brief returns the number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const
    {
      return impl().numBoundarySegments();
    }

    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const
    {
      return impl().size(level,codim);
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return impl().size(codim);
    }

    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return impl().size(level, type);
    }

    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return impl().size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return impl().globalIdSet();
    }

    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const
    {
      return impl().localIdSet();
    }

    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      return impl().levelIndexSet(level);
    }

    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return impl().leafIndexSet();
    }


    /** \brief Create Entity from EntitySeed */
    template<class EntitySeed>
    typename Traits::template Codim<EntitySeed::codimension>::Entity entity(const EntitySeed& seed) const
    {
      return impl().entity(Codim<EntitySeed::codimension>{}, seed);
    }


    /** @name Grid Refinement Methods */
    /*@{*/


    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount)
    {
      impl().globalRefine(refCount);
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
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e)
    {
      return impl().mark(refCount, e);
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark (const typename Traits::template Codim<0>::Entity & e) const
    {
      return impl().getMark(e);
    }

    /** \brief returns true, if at least one entity is marked for adaption */
    bool preAdapt ()
    {
      return impl().preAdapt();
    }

    //! Triggers the grid refinement process
    bool adapt ()
    {
      return impl().adapt();
    }

    /** \brief Clean up refinement markers */
    void postAdapt ()
    {
      return impl().postAdapt();
    }

    /*@}*/


    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize (int codim) const
    {
      return this->leafGridView().overlapSize(codim);
    }

    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize (int codim) const
    {
      return this->leafGridView().ghostSize(codim);
    }

    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize (int level, int codim) const
    {
      return this->levelGridView(level).overlapSize(codim);
    }

    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize (int level, int codim) const
    {
      return this->levelGridView(level).ghostSize(codim);
    }


#if 0
    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    void loadBalance (int strategy, int minlevel, int depth, int maxlevel, int minelement)
    {
      DUNE_THROW(NotImplemented, "Polymorphic::loadBalance()");
    }
#endif


    //! Returns the collective communication object
    const PolymorphicCommunication& comm () const
    {
      return comm_;
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data,
                      InterfaceType iftype,
                      CommunicationDirection dir, int level) const
    {
      PolymorphicCommDataHandle<DataType,Self> dh{data};
      impl().communicate(dh, iftype, dir, level);
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data,
                      InterfaceType iftype,
                      CommunicationDirection dir) const
    {
      PolymorphicCommDataHandle<DataType,Self> dh{data};
      impl().communicate(dh, iftype, dir);
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

    PolymorphicCommunication comm_{};
  };


  namespace Capabilities
  {
    /** \brief has entities for all codimensions
     * \ingroup Polymorphic
     */
    template<int dimension, int dimensionworld, typename ct, int codim>
    struct hasEntity<PolymorphicGrid<dimension, dimensionworld, ct>, codim>
    {
      static const bool v = true; // we assume this
    };

    template<int dimension, int dimensionworld, typename ct, int codim>
    struct hasEntityIterator<PolymorphicGrid<dimension, dimensionworld, ct>, codim>
    {
      static const bool v = true; // we assume this
    };

    /** \brief Polymorphic can communicate
     *  \ingroup Polymorphic
     */
    template<int dimension, int dimensionworld, typename ct, int codim>
    struct canCommunicate<PolymorphicGrid<dimension, dimensionworld, ct>, codim>
    {
      static const bool v = true; // we assume this
    };

    /** \brief has conforming level grids
     * \ingroup Polymorphic
     */
    template<int dimension, int dimensionworld, typename ct>
    struct isLevelwiseConforming<PolymorphicGrid<dimension, dimensionworld, ct>>
    {
      static const bool v = true; // we assume this
    };

    /** \brief has conforming leaf grids
     * \ingroup Polymorphic
     */
    template<int dimension, int dimensionworld, typename ct>
    struct isLeafwiseConforming<PolymorphicGrid<dimension, dimensionworld, ct>>
    {
      static const bool v = true; // we assume this
    };
  } // end namespace Capabilities

} // namespace Dune

#include "io/file/dgfparser/dgfpolymorphic.hh"

#endif // DUNE_GRID_VIRTUALIZEDGRID_HH
