// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRVIRTUALIZED_HH
#define DUNE_VIRTUALIZEDGRVIRTUALIZED_HH

/** \file
 * \brief The VirtualizedGridEntity class
 */

#include <dune/grid/common/grid.hh>
#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

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


  template<int codim, int dim, class GridImp>
  struct VirtualizedGridEntityDefinition
  {
    using Geometry = typename GridImp::template Codim<codim>::Geometry;
    using EntitySeed = typename GridImp::template Codim<codim>::EntitySeed;

    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool equals (const VirtualizedGridEntity<codim, dim, GridImp>& other) const = 0;
      virtual EntitySeed seed () const = 0;
      virtual int level () const = 0;
      virtual PartitionType partitionType () const = 0;
      virtual unsigned int subEntities (unsigned int cc) const = 0;
      virtual Geometry geometry () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      using Wrapped = typename Wrapper::Wrapped;

      bool equals (const VirtualizedGridEntity<codim, dim, GridImp>& other) const final
      {
        return this->get() == Polymorphic::asWrapped<Wrapped>(other);
      }

      EntitySeed seed () const final
      {
        return EntitySeed( this->get().seed() );
      }

      int level () const final
      {
        return this->get().level();
      }

      PartitionType partitionType () const final
      {
        return this->get().partitionType();
      }

      unsigned int subEntities (unsigned int cc) const final
      {
        return this->get().subEntities(cc);
      }

      Geometry geometry () const final
      {
        return Geometry( this->get().geometry() );
      }
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };


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
    public VirtualizedGridEntityDefinition<codim,dim,GridImp>::Base,
    public EntityDefaultImplementation <codim, dim, GridImp, VirtualizedGridEntity>
  {
    using Definition = VirtualizedGridEntityDefinition<codim,dim,GridImp>;
    using Base = typename Definition::Base;

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

  public:

    VirtualizedGridEntity() = default;

    template <class Impl, disableCopyMove<VirtualizedGridEntity,Impl> = 0>
    VirtualizedGridEntity (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    bool equals (const VirtualizedGridEntity& other) const
    {
      return this->asInterface().equals(other);
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return this->asInterface().seed();
    }

    //! level of this element
    int level () const
    {
      return this->asInterface().level();
    }

    /** \brief The partition type for parallel computing
     */
    PartitionType partitionType () const
    {
      return this->asInterface().partitionType();
    }

    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int cc) const
    {
      return this->asInterface().subEntities(cc);
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      return Geometry( this->asInterface().geometry() );
    }
  };



  template<int dim, class GridImp>
  struct VirtualizedGridEntityDefinition<0,dim,GridImp>
  {
    using Geometry = typename GridImp::template Codim<0>::Geometry;
    using LocalGeometry = typename GridImp::template Codim<0>::LocalGeometry;
    using EntitySeed = typename GridImp::template Codim<0>::EntitySeed;
    using LevelIntersectionIterator = VirtualizedGridLevelIntersectionIterator<const GridImp>;
    using LeafIntersectionIterator = VirtualizedGridLeafIntersectionIterator<const GridImp>;
    using HierarchicIterator = VirtualizedGridHierarchicIterator<const GridImp>;

    struct Interface
    {
      virtual ~Interface () = default;
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
      virtual bool isLeaf() const = 0;
      virtual bool hasBoundaryIntersections () const = 0;
      virtual typename GridImp::template Codim<0>::Entity father () const = 0;
      virtual LocalGeometry geometryInFather () const = 0;
      virtual HierarchicIterator hbegin (int maxLevel) const = 0;
      virtual HierarchicIterator hend (int maxLevel) const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      using Wrapped = typename Wrapper::Wrapped;

      bool equals (const VirtualizedGridEntity<0, dim, GridImp>& other) const final
      {
        return this->get() == Polymorphic::asWrapped<Wrapped>(other);
      }

      bool hasFather () const final
      {
        return this->get().hasFather();
      }

      EntitySeed seed () const final
      {
        return EntitySeed( this->get().seed() );
      }

      int level () const final
      {
        return this->get().level();
      }

      PartitionType partitionType () const final
      {
        return this->get().partitionType();
      }

      Geometry geometry () const final
      {
        return Geometry( this->get().geometry() );
      }

      unsigned int subEntities (unsigned int cc) const final
      {
        return this->get().subEntities(cc);
      }

      typename GridImp::template Codim<0>::Entity subEntity0 (int i) const final
      {
        return VirtualizedGridEntity<0, dim, GridImp>( this->get().template subEntity<0>(i) );
      }

      typename GridImp::template Codim<1>::Entity subEntity1 (int i) const final
      {
        return VirtualizedGridEntity<1, dim, GridImp>( this->get().template subEntity<1>(i) );
      }

      typename GridImp::template Codim<dim-1>::Entity subEntityDimMinus1 (int i) const final
      {
        return VirtualizedGridEntity<dim-1, dim, GridImp>( this->get().template subEntity<dim-1>(i) );
      }

      typename GridImp::template Codim<dim>::Entity subEntityDim (int i) const final
      {
        return VirtualizedGridEntity<dim, dim, GridImp>( this->get().template subEntity<dim>(i) );
      }

      bool isLeaf () const final
      {
        return this->get().isLeaf();
      }

      bool hasBoundaryIntersections () const final
      {
        return this->get().hasBoundaryIntersections();
      }

      typename GridImp::template Codim<0>::Entity father () const final
      {
        return VirtualizedGridEntity<0, dim, GridImp>(this->get().father());
      }

      LocalGeometry geometryInFather () const final
      {
        return LocalGeometry(this->get().geometryInFather());
      }

      HierarchicIterator hbegin (int maxLevel) const final
      {
        return VirtualizedGridHierarchicIterator<const GridImp>(this->get().hbegin(maxLevel));
      }

      HierarchicIterator hend (int maxLevel) const final
      {
        return VirtualizedGridHierarchicIterator<const GridImp>(this->get().hend(maxLevel));
      }
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
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
    public VirtualizedGridEntityDefinition<0,dim,GridImp>::Base,
    public EntityDefaultImplementation<0, dim, GridImp, VirtualizedGridEntity>
  {
    using Definition = VirtualizedGridEntityDefinition<0,dim,GridImp>;
    using Base = typename Definition::Base;

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

  public:

    VirtualizedGridEntity() = default;

    template <class Impl, Dune::disableCopyMove<VirtualizedGridEntity,Impl> = 0>
    VirtualizedGridEntity (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {
      // static_assert(Concept::models<typename Definition::Concept,Impl>(),
      //   "Implementation does not model the GridEntity concept.");
    }

    bool equals (const VirtualizedGridEntity& other) const
    {
      return this->asInterface().equals(other);
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return this->asInterface().hasFather();
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return this->asInterface().seed();
    }

    //! Level of this element
    int level () const
    {
      return this->asInterface().level();
    }


    /** \brief The partition type for parallel computing */
    PartitionType partitionType () const
    {
      return this->asInterface().partitionType();
    }


    //! Geometry of this entity
    Geometry geometry () const
    {
      return Geometry( this->asInterface().geometry() );
    }


    /** \brief Return the number of subEntities of codimension codim.
     */
    unsigned int subEntities (unsigned int codim) const
    {
      return this->asInterface().subEntities(codim);
    }


    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template< int cc >
    typename GridImp::template Codim<cc>::Entity subEntity (int i) const
    {
      if constexpr (cc == 0)
        return this->asInterface().subEntity0(i);
      if constexpr (cc == 1)
        return this->asInterface().subEntity1(i);
      if constexpr (cc == dim-1)
        return this->asInterface().subEntityDimMinus1(i);
      if constexpr (cc == dim)
        return this->asInterface().subEntityDim(i);
    }


    //! returns true if Entity has NO children
    bool isLeaf() const
    {
      return this->asInterface().isLeaf();
    }

    // returns true if entity has boundary intersections
    bool hasBoundaryIntersections () const
    {
      return this->asInterface().hasBoundaryIntersections();
    }


    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    typename GridImp::template Codim<0>::Entity father () const
    {
      return this->asInterface().father();
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
      return LocalGeometry( this->asInterface().geometryInFather() );
    }


    /** \brief Inter-level access to son elements on higher levels<=maxlevel.
     * This is provided for sparsely stored nested unstructured meshes.
     * Returns iterator to first son.
     */
    VirtualizedGridHierarchicIterator<const GridImp> hbegin (int maxLevel) const
    {
      return VirtualizedGridHierarchicIterator<const GridImp>( this->asInterface().hbegin(maxLevel) );
    }


    //! Returns iterator to one past the last son
    VirtualizedGridHierarchicIterator<const GridImp> hend (int maxLevel) const
    {
      return VirtualizedGridHierarchicIterator<const GridImp>( this->asInterface().hend(maxLevel) );
    }

  }; // end of VirtualizedGridEntity codim = 0


} // namespace Dune


#endif
