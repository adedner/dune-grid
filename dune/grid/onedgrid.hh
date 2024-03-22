// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_HH
#define DUNE_ONE_D_GRID_HH

#include <tuple>
#include <vector>
#include <list>

#include <dune/common/parallel/communication.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/geometry/type.hh>

/** \file
 * \brief The OneDGrid class
 */

#include "onedgrid/onedgridlist.hh"
#include "onedgrid/nulliteratorfactory.hh"
#include "onedgrid/onedgridentity.hh"
#include "onedgrid/onedgridentityseed.hh"
#include "onedgrid/onedgridintersections.hh"
#include "onedgrid/onedgridintersectioniterators.hh"
#include "onedgrid/onedgridleafiterator.hh"
#include "onedgrid/onedgridviews.hh"
#include "onedgrid/onedgridleveliterator.hh"
#include "onedgrid/onedgridhieriterator.hh"
#include "onedgrid/onedgridindexsets.hh"

namespace Dune {

  template <int dimw, class ct>
  class OneDEmbeddedGrid;

  using OneDGrid = OneDEmbeddedGrid<1,double>;

  /** \brief The type used to for OneDGrid geometries

    If you ever want OneDGrid to use a different type for coordinates,
    you need to change the first argument of AxisAlignedCubeGeometry here.
  */

  template <int dimw, class ctype>
  struct OneDGridFamily
  {
    template <int mydim, int coorddim, class GridImp>
    using OneDGridGeometry = AffineGeometry<ctype, mydim, coorddim>;

    template <int mydim, int coorddim, class GridImp>
    using OneDGridLocalGeometry = AxisAlignedCubeGeometry<ctype, mydim, coorddim>;

    typedef OneDEmbeddedGrid<dimw,ctype> Grid;
    typedef GridTraits<1,     // Grid dimension
                       dimw,  // Dimension of the physical space
        Grid,
        OneDGridGeometry,
        OneDGridEntity,
        OneDGridLevelIterator,
        OneDGridLeafIntersection,
        OneDGridLevelIntersection,
        OneDGridLeafIntersectionIterator,
        OneDGridLevelIntersectionIterator,
        OneDGridHierarchicIterator,
        OneDGridLeafIterator,
        OneDGridLevelIndexSet<const Grid>,
        OneDGridLeafIndexSet<const Grid>,
        OneDGridIdSet<const Grid>,
        unsigned int,
        OneDGridIdSet<const Grid>,
        unsigned int,
        Communication<No_Comm>,
        OneDGridLevelGridViewTraits,
        OneDGridLeafGridViewTraits,
        OneDGridEntitySeed,
        OneDGridLocalGeometry,
        unsigned int,
        std::array<GeometryType,1> >
    Traits;
  };

  //**********************************************************************
  //
  // --OneDGrid
  //
  //**********************************************************************

  /**
     \brief One-dimensional adaptive grid

     [<em> provides \ref Dune::Grid </em>]
     \ingroup GridImplementations
     \ingroup OneDGrid

     This implementation of the grid interface provides one-dimensional
     grids only. The OneDEmbeddedGrid can be nonuniform
     and provides local mesh refinement and coarsening.
   */
  template <int dimw = 1, class ct = double>
  class OneDEmbeddedGrid : public GridDefaultImplementation <1, dimw, ct, OneDGridFamily<dimw,ct>>
  {
    // Grid and world dimension are hardwired in this grid
    constexpr static int dim = 1;
    constexpr static int dimworld = dimw;

    template <int , PartitionIteratorType, class >
    friend class OneDGridLevelIterator;

    friend class OneDGridHierarchicIterator<const OneDEmbeddedGrid>;

    template <int codim_, int dim_, class GridImp_>
    friend class OneDGridEntity;
    friend class OneDGridHierarchicIterator<OneDEmbeddedGrid>;
    friend class OneDGridLeafIntersection<const OneDEmbeddedGrid>;
    friend class OneDGridLevelIntersection<const OneDEmbeddedGrid>;
    friend class OneDGridLeafIntersectionIterator<const OneDEmbeddedGrid>;
    friend class OneDGridLevelIntersectionIterator<const OneDEmbeddedGrid>;

    friend class OneDGridLevelIndexSet<const OneDEmbeddedGrid>;
    friend class OneDGridLeafIndexSet<const OneDEmbeddedGrid>;
    friend class OneDGridIdSet<const OneDEmbeddedGrid>;

    template <int codim_, PartitionIteratorType PiType_, class GridImp_>
    friend class OneDGridLeafIterator;

    friend class OneDGridLeafGridView<const OneDEmbeddedGrid>;
    friend class OneDGridLevelGridView<const OneDEmbeddedGrid>;

    template <class GridType_>
    friend class GridFactory;

    template<int codim_, int dim_, class GridImp_, template<int,int,class> class EntityImp_>
    friend class Entity;

    /** \brief Default constructor for the GridFactory */
    OneDEmbeddedGrid();

    // **********************************************************
    // The Interface Methods
    // **********************************************************

  public:

    /** \brief The type used to store coordinates
     */
    typedef ct ctype;

    /** \brief GridFamily of OneDGrid */
    typedef OneDGridFamily<dimw,ct> GridFamily;

    //Provides the standard grid types
    typedef typename GridFamily::Traits Traits;

    typedef typename Traits::template Codim<0>::Entity::Geometry::GlobalCoordinate GlobalCoordinate;

    /** \brief Constructor with an explicit set of coordinates */
    template <class C>
    OneDEmbeddedGrid(const std::vector<C>& coords);

    /** \brief Constructor for a uniform grid */
    template <int d = dimw, std::enable_if_t<(d == 1), int> = 0>
    OneDEmbeddedGrid(int numElements, const ctype& leftBoundary, const ctype& rightBoundary)
      : OneDEmbeddedGrid(numElements, GlobalCoordinate(leftBoundary), GlobalCoordinate(rightBoundary))
    {}

    /** \brief Constructor for a uniform grid */
    OneDEmbeddedGrid(int numElements, const GlobalCoordinate& leftBoundary, const GlobalCoordinate& rightBoundary);

    //! Destructor
    ~OneDEmbeddedGrid();

    /** \brief Return maximum level defined in this grid.

       Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {return entityImps_.size()-1;}

    /** \brief Create an Entity from an EntitySeed */
    template <typename Seed>
    static typename Traits::template Codim<Seed::codimension>::Entity
    entity(const Seed& seed)
    {
      const int codim = Seed::codimension;
      return typename Traits::template Codim<codim>::Entity(OneDGridEntity<codim,dim,const OneDEmbeddedGrid>(seed.impl().target()));
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      switch (codim)
      {
        case 0:
          return elements(level).size();
        case 1:
          return vertices(level).size();
        default:
          return 0;
      }
    }



    //! number of leaf entities per codim in this process
    int size (int codim) const
    {
      return leafIndexSet().size(codim);
    }

    //! number of entities per level and geometry type in this process
    int size (int level, GeometryType type) const
    {
      // There is only one type for each codim
      return size(level,1-type.dim());
    }

    //! number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }

    /** \brief Return the number of coarse grid boundary segments.

       For this grid implementation, the return value is always 2, because only connected domains
       are supported, and then the coarse grid boundary consists of two points.
     */
    size_t numBoundarySegments() const
    {
      return 2;
    }

    /** \brief Get the set of global ids */
    const typename Traits::GlobalIdSet& globalIdSet() const
    {
      return idSet_;
    }

    /** \brief Get the set of local ids */
    const typename Traits::LocalIdSet& localIdSet() const
    {
      return idSet_;
    }

    /** \brief Get an index set for the given level */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (! levelIndexSets_[level]) {
        levelIndexSets_[level] =
          new OneDGridLevelIndexSet<const OneDEmbeddedGrid>(*this, level);
        levelIndexSets_[level]->update();
      }

      return * levelIndexSets_[level];
    }

    /** \brief Get an index set for the leaf level */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }


    /** \brief Mark entity for refinement
     *
     * \param refCount if >0 mark for refinement, if <0 mark for coarsening
     * \param e Entity to the entity you want to mark
     *
     * \return True, if marking was successful
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity& e );

    /** \brief return current adaptation marker of given entity

        \param e Entity to the entity you want to mark

        \return int current adaptation marker of entity pointer e
     */
    int getMark(const typename Traits::template Codim<0>::Entity& e ) const;

    //! Does nothing except return true if some element has been marked for refinement
    bool preAdapt();

    //! Triggers the grid refinement process
    bool adapt();

    /** \brief Adaptation post-processing: Reset all adaptation state flags */
    void postAdapt();

    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    /** \brief The different forms of grid refinement supported by OneDGrid */
    enum RefinementType {
      /** \brief New level consists only of the refined elements */
      LOCAL,
      /** \brief New level consists of the refined elements and the unrefined ones, too */
      COPY
    };

    /** \brief Sets the type of grid refinement */
    void setRefinementType(RefinementType type) {
      refinementType_ = type;
    }

    /** \brief Does one uniform refinement step
     *
     * \param refCount I don't know what this is good for.  It doesn't
     *        actually do anything.
     */
    void globalRefine(int refCount);

    // dummy parallel functions

    const Communication<No_Comm> &comm () const
    {
      return ccobj;
    }


    /** \brief Communicate data of level gridView */
    template <class DataHandle>
    void communicate (DataHandle& /*handle*/, InterfaceType /*iftype*/,
                      CommunicationDirection /*dir*/, int /*level*/) const
    {
      DUNE_THROW(Dune::NotImplemented, "communicate() for OneDGrid not implemented");
    }

    /** \brief Communicate data of leaf gridView */
    template <class DataHandle>
    void communicate (DataHandle& /*handle*/, InterfaceType /*iftype*/,
                      CommunicationDirection /*dir*/) const
    {
      DUNE_THROW(Dune::NotImplemented, "communicate() for OneDGrid not implemented");
    }


  private:

    /** \brief Get vertex lists directly -- makes the code more readable */
    OneDGridList<OneDEntityImp<0,dimw,ct> >& vertices(int level) {
      return std::get<0>(entityImps_[level]);
    }

    /** \brief Get vertex lists directly -- makes the code more readable */
    const OneDGridList<OneDEntityImp<0,dimw,ct> >& vertices(int level) const {
      return std::get<0>(entityImps_[level]);
    }

    /** \brief Get element lists directly -- makes the code more readable */
    OneDGridList<OneDEntityImp<1,dimw,ct> >& elements(int level) {
      return std::get<1>(entityImps_[level]);
    }

    /** \brief Get element lists directly -- makes the code more readable */
    const OneDGridList<OneDEntityImp<1,dimw,ct> >& elements(int level) const {
      return std::get<1>(entityImps_[level]);
    }

    Communication<No_Comm> ccobj;

    /** \brief Update all indices and ids */
    void setIndices();

    /** \brief Produce an entity id that has not been used in this grid before.
     */
    unsigned int getNextFreeId()
    {
      return freeIdCounter_++;
    }

    //! The type of grid refinement currently in use
    RefinementType refinementType_;

    typename OneDGridList<OneDEntityImp<0,dimw,ct> >::iterator getLeftUpperVertex(const OneDEntityImp<1,dimw,ct>* eIt);

    typename OneDGridList<OneDEntityImp<0,dimw,ct> >::iterator getRightUpperVertex(const OneDEntityImp<1,dimw,ct>* eIt);

    /** \brief Returns an iterator to the first element on the left of
        the input element which has sons.
     */
    typename OneDGridList<OneDEntityImp<1,dimw,ct> >::iterator getLeftNeighborWithSon(typename OneDGridList<OneDEntityImp<1,dimw,ct> >::iterator eIt);

    // The vertices and elements of the grid hierarchy
    std::vector<std::tuple<OneDGridList<OneDEntityImp<0,dimw,ct> >,
            OneDGridList<OneDEntityImp<1,dimw,ct> > > > entityImps_;

    // Our set of level indices
    mutable std::vector<OneDGridLevelIndexSet<const OneDEmbeddedGrid>* > levelIndexSets_;

    OneDGridLeafIndexSet<const OneDEmbeddedGrid> leafIndexSet_;

    OneDGridIdSet<const OneDEmbeddedGrid> idSet_;

    // Every entity gets a unique id, unless it is a copy of an entity on a coarser level.
    // This is the counter that we use to create the unique id.
    unsigned int freeIdCounter_;

    /** Since a OneDGrid is one-dimensional and connected, there can only be two possible numberings
        of the boundary segments.  Either the left one is '0' and the right one is '1' or the reverse.
        This flag stores which is the case. */
    bool reversedBoundarySegmentNumbering_;

  }; // end Class OneDEmbeddedGrid

  extern template class OneDEmbeddedGrid<1,double>;
  extern template OneDEmbeddedGrid<1,double>::OneDEmbeddedGrid(const std::vector<double>&);
  extern template OneDEmbeddedGrid<1,double>::OneDEmbeddedGrid(const std::vector<FieldVector<double,1>>&);
  extern template class OneDEmbeddedGrid<2,double>;
  extern template OneDEmbeddedGrid<2,double>::OneDEmbeddedGrid(const std::vector<FieldVector<double,2>>&);
  extern template class OneDEmbeddedGrid<3,double>;
  extern template OneDEmbeddedGrid<3,double>::OneDEmbeddedGrid(const std::vector<FieldVector<double,3>>&);

  namespace Capabilities
  {
    /** \struct hasBackupRestoreFacilities
       \ingroup OneDGrid
     */

    /** \struct IsUnstructured
       \ingroup OneDGrid
     */

    /** \brief OneDGrid has only one geometry type for codim 0 entities
       \ingroup OneDGrid
     */
    template<int dimw, class ct>
    struct hasSingleGeometryType< OneDEmbeddedGrid<dimw,ct> >
    {
      static const bool v = true;
      static const unsigned int topologyId = GeometryTypes::cube(1).id();
    };


    /** \brief OneDGrid has entities for all codimension
       \ingroup OneDGrid
     */
    template<int dimw, class ct, int cdim>
    struct hasEntity< OneDEmbeddedGrid<dimw,ct>, cdim >
    {
      static const bool v = true;
    };

    /**
     * \brief OneDGrid can iterate over all codimensions
     * \ingroup OneDGrid
     **/
    template<int dimw, class ct, int codim>
    struct hasEntityIterator<OneDEmbeddedGrid<dimw,ct>, codim>
    {
      static const bool v = true;
    };

    /** \brief OneDGrid is levelwise conforming
       \ingroup OneDGrid
     */
    template<int dimw, class ct>
    struct isLevelwiseConforming< OneDEmbeddedGrid<dimw,ct> >
    {
      static const bool v = true;
    };

    /** \brief OneDGrid is leafwise conforming
       \ingroup OneDGrid
     */
    template<int dimw, class ct>
    struct isLeafwiseConforming< OneDEmbeddedGrid<dimw,ct> >
    {
      static const bool v = true;
    };

  }

} // namespace Dune

// Include the GridFactory specialization for OneDGrid, so everybody
// who includes the grid also gets the factory.  Since OneDGrid is
// not a template class, it needs to be a complete type before
// GridFactory<OneDGrid> can be defined.  This is why the #include-
// directive is at _the end_ of this file.
#include <dune/grid/onedgrid/onedgridfactory.hh>


#endif
