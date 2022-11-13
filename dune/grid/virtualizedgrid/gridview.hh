// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZED_GRIDVIEW_HH
#define DUNE_VIRTUALIZED_GRIDVIEW_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/grid/virtualizedgrid/indexsets.hh>
#include <dune/grid/virtualizedgrid/intersections.hh>
#include <dune/grid/virtualizedgrid/intersectioniterator.hh>

namespace Dune
{

  //! Forward declarations
  template< class GridImp >
  class VirtualizedLevelGridView;

  template< class GridImp >
  class VirtualizedLeafGridView;


  // VirtualizedGridLevelViewTraits
  // ------------------------------

  template< class GridImp >
  struct VirtualizedGridLevelViewTraits
  {
    typedef VirtualizedLevelGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid::Traits::LevelIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid::Traits::LevelIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid::Traits::LevelIntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid::Traits::Communication Communication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< cd >::template Partition< All_Partition >::LevelIterator Iterator;

      typedef typename Grid::Traits::template Codim< cd >::Entity Entity;

      typedef typename Grid::template Codim< cd >::Geometry Geometry;
      typedef typename Grid::template Codim< cd >::LocalGeometry LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid::template Codim< cd >::template Partition< pit >::LevelIterator Iterator;
      };
    };

    constexpr static bool conforming = Capabilities::isLevelwiseConforming< Grid >::v;
  };

  // VirtualizedGridView
  // -------------------

  template< class GridImp >
  class VirtualizedLevelGridView
  {
    typedef VirtualizedLevelGridView< GridImp > ThisType;

  public:
    typedef VirtualizedGridLevelViewTraits< GridImp > Traits;

    typedef typename Traits::Grid Grid;

    typedef typename Traits::IndexSet IndexSet;

    typedef typename Traits::Intersection Intersection;

    typedef typename Traits::IntersectionIterator IntersectionIterator;

    typedef typename Traits::Communication Communication;

    template< int codim >
    struct Codim : public Traits::template Codim< codim > {};

    static constexpr bool conforming = Traits::conforming;

    VirtualizedLevelGridView ( const Grid &grid, int level )
      : grid_( &grid ), level_( level )
    {}

    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    const IndexSet &indexSet () const
    {
      return grid().levelIndexSet( level_ );
    }

    bool isConforming() const { return bool(Traits::conforming); }

    int size ( int codim ) const
    {
      return grid().size( level_, codim );
    }

    int size ( const GeometryType &type ) const
    {
      return grid().size( level_, type );
    }

    template< int codim >
    typename Codim< codim >::Iterator begin () const
    {
      return grid().template lbegin< codim, All_Partition >( level_ );
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator begin () const
    {
      return grid().template lbegin< codim, pit >( level_ );
    }

    template< int codim >
    typename Codim< codim >::Iterator end () const
    {
      return grid().template lend< codim, All_Partition >( level_ );
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator end () const
    {
      return grid().template lend< codim, pit >( level_ );
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      return grid().ilevelbegin( entity );
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      return grid().ilevelend( entity );
    }

    const Communication &comm () const
    {
      return grid().comm();
    }

    int overlapSize ( int codim ) const
    {
      return grid().overlapSize(level_, codim );
    }

    int ghostSize ( int codim ) const
    {
      return grid().ghostSize(level_, codim );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
      return grid().communicate( dataHandle, interface, direction, level_ );
    }

  private:
    const Grid *grid_;
    int level_;
  };


  // VirtualizedGridLeafViewTraits
  // -----------------------------

  template< class GridImp >
  struct VirtualizedGridLeafViewTraits
  {
    typedef VirtualizedLeafGridView< GridImp > GridViewImp;

    /** \brief type of the grid */
    typedef typename std::remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid::Traits::LeafIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid::Traits::LeafIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid::Traits::LeafIntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid::Traits::Communication Communication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< cd >::template Partition< All_Partition >::LeafIterator Iterator;

      typedef typename Grid::Traits::template Codim< cd >::Entity Entity;

      typedef typename Grid::template Codim< cd >::Geometry Geometry;
      typedef typename Grid::template Codim< cd >::LocalGeometry LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template <PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid::template Codim< cd >::template Partition< pit >::LeafIterator Iterator;
      };
    };

    constexpr static bool conforming = Capabilities::isLeafwiseConforming< Grid >::v;
  };


  // VirtualizedGridView
  // -------------------

  template< class GridImp >
  class VirtualizedLeafGridView
  {
    typedef VirtualizedLeafGridView< GridImp > ThisType;

  public:
    typedef VirtualizedGridLeafViewTraits< GridImp > Traits;

    typedef typename Traits::Grid Grid;

    typedef typename Traits::IndexSet IndexSet;

    typedef typename Traits::Intersection Intersection;

    typedef typename Traits::IntersectionIterator IntersectionIterator;

    typedef typename Traits::Communication Communication;

    template< int codim >
    struct Codim : public Traits::template Codim< codim > {};

    static constexpr bool conforming = Traits::conforming;

    VirtualizedLeafGridView ( const Grid &grid )
    : grid_( &grid )
    {}

    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    const IndexSet &indexSet () const
    {
      return grid().leafIndexSet();
    }

    bool isConforming() const { return bool(Traits::conforming); }

    int size ( int codim ) const
    {
      return grid().size( codim );
    }

    int size ( const GeometryType &type ) const
    {
      return grid().size( type );
    }

    template< int codim >
    typename Codim< codim >::Iterator begin () const
    {
      return grid().template leafbegin< codim, All_Partition >();
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator begin () const
    {
      return grid().template leafbegin< codim, pit >();
    }

    template< int codim >
    typename Codim< codim >::Iterator end () const
    {
      return grid().template leafend< codim, All_Partition >();
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator end () const
    {
      return grid().template leafend< codim, pit >();
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      return grid().ileafbegin( entity );
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      return grid().ileafend( entity );
    }

    const Communication &comm () const
    {
      return grid().comm();
    }

    int overlapSize ( int codim ) const
    {
      return grid().overlapSize( codim );
    }

    int ghostSize ( int codim ) const
    {
      return grid().ghostSize( codim );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                      InterfaceType interface,
                      CommunicationDirection direction ) const
    {
      return grid().communicate( dataHandle, interface, direction );
    }

  private:
    const Grid *grid_;
  };

} // namespace Dune

#endif // #ifndef DUNE_VIRTUALIZEDGRID_GRIDVIEW_HH
