// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_VIRTUALIZEDGRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_HH

/** \file
 * \brief The VirtualizedGrid class
 */

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// The components of the VirtualizedGrid interface
#include "virtualizedgrid/grid.hh"

namespace Dune
{
  // Forward declaration
  template<int dimension, int dimensionworld, typename ct>
  class VirtualizedGrid;

  template<int dimension, int dimensionworld, typename ct>
  struct VirtualizedGridFamily
  {
    struct Traits
    {
      /** \brief The type that implements the grid. */
      typedef VirtualizedGrid< dimension, dimensionworld, ct > Grid;

      using IntersectionImp = VirtualizedGridIntersection< const Grid >;
      using IntersectionIteratorImp = VirtualizedGridIntersectionIterator<const Grid>;

      /** \brief The type of the intersection at the leafs of the grid. */
      using LeafIntersection = Dune::Intersection<const Grid, IntersectionImp>;

      /** \brief The type of the intersection at the levels of the grid. */
      using LevelIntersection = Dune::Intersection<const Grid, IntersectionImp>;

      /** \brief The type of the intersection iterator at the leafs of the grid. */
      using LeafIntersectionIterator = Dune::IntersectionIterator<const Grid, IntersectionIteratorImp, IntersectionImp>;

      /** \brief The type of the intersection iterator at the levels of the grid. */
      using LevelIntersectionIterator = Dune::IntersectionIterator<const Grid, IntersectionIteratorImp, IntersectionImp>;

      /** \brief The type of the  hierarchic iterator. */
      typedef Dune::EntityIterator< 0, const Grid, VirtualizedGridEntityIterator<0,const Grid> > HierarchicIterator;

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
          using LevelIterator = Dune::EntityIterator<cd, const Grid, VirtualizedGridEntityIterator<cd, const Grid> >;
          /** \brief The type of the iterator over the leaf entities of this codim on this partition. */
          using LeafIterator = Dune::EntityIterator< cd, const Grid, VirtualizedGridEntityIterator<cd, const Grid> >;
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
      typedef IndexSet< const Grid, VirtualizedGridIndexSet< const Grid > > LevelIndexSet;
      /** \brief The type of the leaf index set. */
      typedef IndexSet< const Grid, VirtualizedGridIndexSet< const Grid > > LeafIndexSet;
      /** \brief The type of the global id set. */
      typedef IdSet< const Grid, VirtualizedGridIdSet< const Grid >, VirtualizedGridIdType> GlobalIdSet;
      /** \brief The type of the local id set. */
      typedef IdSet< const Grid, VirtualizedGridIdSet< const Grid >, VirtualizedGridIdType> LocalIdSet;

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

  template<int dimension, int dimensionworld, typename ct>
  class VirtualizedGrid
    : public VirtualizedGridDefinition<VirtualizedGrid<dimension,dimensionworld,ct>,dimension,dimensionworld,ct>::Base
    , public GridDefaultImplementation<dimension, dimensionworld, ct, VirtualizedGridFamily<dimension, dimensionworld, ct>>
  {
    using Self = VirtualizedGrid<dimension, dimensionworld, ct>;
    using Definition = VirtualizedGridDefinition<Self,dimension,dimensionworld,ct>;
    using Base = typename Definition::Base;

    template <PartitionIteratorType p>
    using _Partition = typename Definition::template _Partition<p>;

  public:
    //! type of the used GridFamily for this grid
    typedef VirtualizedGridFamily<dimension, dimensionworld, ct> GridFamily;

    //! the Traits
    typedef typename GridFamily::Traits Traits;

    //! The type used to store coordinates
    typedef ct ctype;

  public:

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    /** \brief Constructor
     *
     * \param grid The grid hold by the VirtualizedGrid
     */
    template <class Impl, disableCopyMove<VirtualizedGrid,Impl> = 0>
    VirtualizedGrid (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    /** \brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const
    {
      return this->asInterface().maxLevel();
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const
    {
      return this->asInterface().lbegin(Codim<codim>{},level);
    }

    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const
    {
      return this->asInterface().lend(Codim<codim>{},level);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lbegin (int level) const
    {
      return this->asInterface().lbegin(Codim<codim>{},_Partition<pitype>{},level);
    }

    //! one past the end on this level
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator lend (int level) const
    {
      return this->asInterface().lend(Codim<codim>{},_Partition<pitype>{},level);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin () const
    {
      return this->asInterface().leafbegin(Codim<codim>{});
    }

    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend () const
    {
      return this->asInterface().leafend(Codim<codim>{});
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafbegin () const
    {
      return this->asInterface().leafbegin(Codim<codim>{},_Partition<pitype>{});
    }

    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType pitype>
    typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator leafend () const
    {
      return this->asInterface().leafend(Codim<codim>{},_Partition<pitype>{});
    }


    virtual typename Traits::LevelIntersectionIterator ilevelbegin (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return this->asInterface().ilevelbegin( entity );
    }

    virtual typename Traits::LevelIntersectionIterator ilevelend (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return this->asInterface().ilevelend( entity );
    }

    virtual typename Traits::LeafIntersectionIterator ileafbegin (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return this->asInterface().ileafbegin( entity );
    }

    virtual typename Traits::LeafIntersectionIterator ileafend (const typename Traits::template Codim<0>::Entity& entity) const
    {
      return this->asInterface().ileafend( entity );
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const
    {
      return this->asInterface().size(level,codim);
    }

    /** \brief returns the number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const
    {
      return this->asInterface().numBoundarySegments();
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return this->asInterface().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const
    {
      return this->asInterface().size(level, type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return this->asInterface().size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return this->asInterface().globalIdSet();
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const
    {
      return this->asInterface().localIdSet();
    }


    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      return this->asInterface().levelIndexSet(level);
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return this->asInterface().leafIndexSet();
    }


    /** \brief Create Entity from EntitySeed */
    template<class EntitySeed>
    typename Traits::template Codim<EntitySeed::codimension>::Entity entity(const EntitySeed& seed) const
    {
      return this->asInterface().entity(Codim<EntitySeed::codimension>{}, seed);
    }


    /** @name Grid Refinement Methods */
    /*@{*/


    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount)
    {
      this->asInterface().globalRefine(refCount);
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
      return this->asInterface().mark(refCount, e);
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark (const typename Traits::template Codim<0>::Entity & e) const
    {
      return this->asInterface().getMark(e);
    }

    /** \brief returns true, if at least one entity is marked for adaption */
    bool preAdapt ()
    {
      return this->asInterface().preAdapt();
    }


    //! Triggers the grid refinement process
    bool adapt ()
    {
      return this->asInterface().adapt();
    }

    /** \brief Clean up refinement markers */
    void postAdapt ()
    {
      return this->asInterface().postAdapt();
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
      DUNE_THROW(NotImplemented, "VirtualizedGrid::loadBalance()");
    }
#endif


    //! Returns the collective communication object
    const VirtualizedCommunication& comm () const
    {
      return ccobj;
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data,
                      InterfaceType iftype,
                      CommunicationDirection dir, int level) const
    {
      VirtualizedCommDataHandle<DataType,Self> dh{data};
      this->asInterface().communicate(dh, iftype, dir, level);
    }

    //! The new communication interface
    template<class DataHandleImp, class DataType>
    void communicate (CommDataHandleIF<DataHandleImp, DataType>& data,
                      InterfaceType iftype,
                      CommunicationDirection dir) const
    {
      VirtualizedCommDataHandle<DataType,Self> dh{data};
      this->asInterface().communicate(dh, iftype, dir);
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

  private:
    VirtualizedCommunication ccobj;
  };


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
