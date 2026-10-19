// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKITERATORS_HH
#define DUNE_GRID_TEST_CHECKITERATORS_HH

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/test/iteratortest.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/test/periodicmap.hh>

template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& ){}
};

// CheckCodimIterators
// -------------------

template< class GridView, int codim,
    bool hasEntity = Dune::Capabilities::hasEntityIterator< typename GridView::Grid, codim >::v >
struct CheckCodimIterators;

template< class GridView, int codim >
struct CheckCodimIterators< GridView, codim, false >
{
  static void apply ( const GridView &gridView,
                      const Dune::GridCheck::PeriodicMap< GridView > & = Dune::GridCheck::noPeriodicMap< GridView >() )
  {
    std::cerr << "Warning: Not checking iterators for codimension " << codim
              << ", because the corresponding entities are not implemented." << std::endl;
  }
};

template< class GridView, int codim >
struct CheckCodimIterators< GridView, codim, true >
{
  static void apply ( const GridView &gridView,
                      const Dune::GridCheck::PeriodicMap< GridView > &periodicMap = Dune::GridCheck::noPeriodicMap< GridView >() );
};


// CheckIterators
// --------------

template< class GridView >
class CheckIterators
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gridView,
                      const Dune::GridCheck::PeriodicMap< GridView > &periodicMap = Dune::GridCheck::noPeriodicMap< GridView >() )
  {
    std::cout << "Checking iterators for higher codimension..." << std::endl;
    Dune::Hybrid::forEach( std::make_index_sequence< GridView::dimension >{}, [ & ]( auto i ){ CheckCodim< i+1 >::apply( gridView, periodicMap ); } );
  }
};


// CheckIterators::CheckCodim
// --------------------------

template< class GridView >
template< int codim >
struct CheckIterators< GridView >::CheckCodim
{
  static void apply ( const GridView &gridView, const Dune::GridCheck::PeriodicMap< GridView > &periodicMap )
  {
    return CheckCodimIterators< GridView, codim >::apply( gridView, periodicMap );
  }
};


// Implementation of CheckCodimIterators
// -------------------------------------

template< class GridView, int codim >
inline void CheckCodimIterators< GridView, codim, true >
::apply ( const GridView &gridView, const Dune::GridCheck::PeriodicMap< GridView > &periodicMap )
{
  typedef typename GridView::Grid::Traits::LocalIdSet::IdType IdType;

  const auto &idSet = gridView.grid().localIdSet();

  // periodic copies of the same entity are different entities of the grid view
  // (with different indices), but they share the same id
  const Dune::GridCheck::PeriodicIdentification< GridView > periodic( gridView, codim, periodicMap );
  using Key = typename Dune::GridCheck::PeriodicIdentification< GridView >::Key;

  std::map< IdType, int > count;
  std::map< Key, int > keyCount;

  const auto codimEnd = gridView.template end< codim >();
  for( auto it = gridView.template begin< codim >(); it != codimEnd; ++it )
  {
    ++count[ idSet.id( *it ) ];
    ++keyCount[ periodic.key( *it ) ];
  }

  const auto elementEnd = gridView.template end< 0 >();
  for( auto it = gridView.template begin< 0 >(); it != elementEnd; ++it )
  {
    const auto &entity = *it;
    for( std::size_t i = 0; i < entity.subEntities(codim); ++i )
    {
      // each entity is visited exactly once ...
      const Key key = periodic.subKey( entity, i, codim );
      if( keyCount[ key ] != 1 )
      {
        std::cerr << "Error: Codim " << codim << " iterator"
                  << " visited entity with index " << key.second
                  << " " << keyCount[ key ] << " times." << std::endl;
        assert( keyCount[ key ] == 1 );
      }

      // ... and so is each id, once for each periodic copy on this process
      IdType id = idSet.subId( entity, i, codim );
      const int copies = periodic.members( periodic.classOf( key ) ).size();
      if( count[ id ] != copies )
      {
        std::cerr << "Error: Codim " << codim << " iterator"
                  << " visited entity " << id
                  << " " << count[ id ] << " times"
                  << " (expected " << copies << " periodic copies)." << std::endl;
        assert( count[ id ] == copies );
      }
    }
  }

  // check forward iterator semantics
  typedef typename GridView::template Codim<codim>::Entity Entity;
  NoopFunctor<Entity> op;
  if(0 != testForwardIterator(gridView.template begin<codim>(),
                               gridView.template end<codim>(), op))
    DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

}



/** \brief check the codim iterators of a grid view
 *
 * \param periodicMap maps each coordinate to its representative in the periodic
 *                    unit cell, see dune/grid/test/periodicmap.hh
 */
template< class GridView >
inline void checkIterators ( const GridView &gridView,
                             const Dune::GridCheck::PeriodicMap< GridView > &periodicMap = Dune::GridCheck::noPeriodicMap< GridView >() )
{
  CheckIterators< GridView >::apply( gridView, periodicMap );
}

#endif // DUNE_GRID_TEST_CHECKITERATORS_HH
