// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_TEST_CHECKPARTITION_HH
#define DUNE_GRID_TEST_CHECKPARTITION_HH

#include <cstddef>

#include <bitset>
#include <set>
#include <map>
#include <tuple>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/entitycommhelper.hh>
#include <dune/grid/test/periodicmap.hh>

template< Dune::PartitionIteratorType pitype >
struct PartitionFilter;

template<>
struct PartitionFilter< Dune::Interior_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType == Dune::InteriorEntity);
  }
};

template<>
struct PartitionFilter< Dune::InteriorBorder_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType == Dune::InteriorEntity) || (partitionType == Dune::BorderEntity);
  }
};

template<>
struct PartitionFilter< Dune::Overlap_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType != Dune::FrontEntity) && (partitionType != Dune::GhostEntity);
  }
};

template<>
struct PartitionFilter< Dune::OverlapFront_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType != Dune::GhostEntity);
  }
};

template<>
struct PartitionFilter< Dune::All_Partition >
{
  static bool contains ([[maybe_unused]] const Dune::PartitionType partitionType)
  {
    return true;
  }
};

template<>
struct PartitionFilter< Dune::Ghost_Partition >
{
  static bool contains ( const Dune::PartitionType partitionType )
  {
    return (partitionType == Dune::GhostEntity);
  }
};



inline bool possibleSubPartitionType ( Dune::PartitionType ept, Dune::PartitionType pt )
{
  switch( ept )
  {
  case Dune::InteriorEntity :
    return (pt == Dune::InteriorEntity) || (pt == Dune::BorderEntity);
  case Dune::OverlapEntity :
    return (pt == Dune::BorderEntity) || (pt == Dune::OverlapEntity) || (pt == Dune::FrontEntity);
  case Dune::GhostEntity :
    return (pt == Dune::BorderEntity) || (pt == Dune::FrontEntity) || (pt == Dune::GhostEntity);
  default :
    std::cerr << "Error: Codimension 0 entity cannot be of partition type " << ept << "." << std::endl;
    return false;
  }
}



template< class GridView, Dune::PartitionIteratorType pitype >
class CheckPartitionType
{
  template< int codim >
  struct CheckCodim;

public:
  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking iterators for " << pitype << "..." << std::endl;
    Dune::Hybrid::forEach( std::make_index_sequence< GridView::dimension+1 >{}, [ & ]( auto i ){ CheckCodim< i >::apply( gridView ); } );
  }
};


template< class GridView, Dune::PartitionIteratorType pitype >
template< int codim >
struct CheckPartitionType< GridView, pitype >::CheckCodim
{
  template< class IdSet >
  static void check ( const std::true_type &, const GridView &gridView,
                      const IdSet &idSet )
  {
    typedef typename GridView::template Codim< codim >::template Partition< pitype >::Iterator Iterator;
    typedef typename GridView::template Codim< 0 >::template Partition< Dune::All_Partition >::Iterator AllIterator;

    // entities are identified by geometry type and index, as periodic copies
    // of an entity are different entities, but share the same id
    const auto &indexSet = gridView.indexSet();
    // (the type is stored as GlobalGeometryTypeIndex, std::pair would compare
    // GeometryTypes by their exact topology id)
    typedef std::pair< std::size_t, typename GridView::IndexSet::IndexType > Key;

    typedef std::map< Key, Dune::PartitionType > Map;
    typedef typename Map::iterator MapIterator;
    Map map;

    try {
      const Iterator end = gridView.template end< codim, pitype >();
      for( Iterator it = gridView.template begin< codim, pitype >(); it != end; ++it )
      {
        Dune::PartitionType pt = it->partitionType();
        if( !PartitionFilter< pitype >::contains( pt ) )
        {
          std::cerr << "Error: Codim " << codim << " iterator for the " << pitype
                    << " visited entity " << idSet.id( *it )
                    << " with partition type " << pt << "." << std::endl;
        }
        map[ Key( Dune::GlobalGeometryTypeIndex::index( it->type() ), indexSet.index( *it ) ) ] = pt;
      }

      const AllIterator allEnd = gridView.template end< 0, Dune::All_Partition >();
      for( AllIterator it = gridView.template begin< 0, Dune::All_Partition >(); it != allEnd; ++it )
      {
        Dune::PartitionType ept = it->partitionType();
        const int subEntities = it->subEntities(codim);
        for( int i = 0; i < subEntities; ++i )
        {
          const auto subEntity = it->template subEntity< codim >( i );
          Dune::PartitionType pt = subEntity.partitionType();
          if( !possibleSubPartitionType( ept, pt ) )
          {
            std::cerr << "Error: Codim " << codim << " entity " << idSet.subId( *it, i, codim )
                      << " with partition type " << pt << " is a subentity of entity "
                      << idSet.id( *it ) << " with partition type " << ept << "." << std::endl;
          }

          const MapIterator mapIt = map.find( Key( Dune::GlobalGeometryTypeIndex::index( subEntity.type() ),
                                                   indexSet.subIndex( *it, i, codim ) ) );
          if( mapIt == map.end() )
          {
            if( PartitionFilter< pitype >::contains( pt ) )
            {
              std::cerr << "Error: Codim " << codim << " entity " << idSet.subId( *it, i, codim )
                        << " with partition type " << pt << " is not visited by codim " << codim
                        << " iterator for the " << pitype << "." << std::endl;
            }
          }
          else
          {
            if( pt != mapIt->second )
            {
              std::cerr << "Error: Codim " << codim << " entity " << idSet.subId( *it, i, codim )
                        << " with partition type " << pt << " reported partition type "
                        << mapIt->second << " when accessed with the codim " << codim
                        << " iterator for the " << pitype << "." << std::endl;
            }
          }
        }
      }
    }
    catch( const Dune::Exception &exception )
    {
      std::cerr << "Error: Caught  exception when testing partition iterator for the " << pitype
                << "(" << exception << ")." << std::endl;
    }
  }

  template< class IdSet >
  static void check (const std::false_type &,
                     [[maybe_unused]] const GridView &gridView,
                     [[maybe_unused]] const IdSet &idSet)
  {}

  static void apply ( const GridView &gridView )
  {
    std::integral_constant<
        bool, Dune::Capabilities::hasEntity< typename GridView::Grid, codim >::v && Dune::Capabilities::hasEntityIterator< typename GridView::Grid, codim >::v
        > capabilityVariable;
    check( capabilityVariable, gridView, gridView.grid().localIdSet() );
  }
};



template< class GridView, Dune::InterfaceType iftype  >
class CheckPartitionDataHandle
  : public Dune::CommDataHandleIF< CheckPartitionDataHandle< GridView, iftype >, int >
{
  static const int dimension = GridView::dimension;

  typedef typename GridView::Grid Grid;
  typedef typename GridView::IndexSet IndexSet;

  // entities are identified by codim, geometry type and index, as periodic
  // copies of an entity are different entities, but share the same id
  typedef std::tuple< int, std::size_t, typename IndexSet::IndexType > Key;
  typedef std::set< Key > CommSet;

  template< class Entity >
  Key key ( const Entity &entity ) const
  {
    return Key( Entity::codimension, Dune::GlobalGeometryTypeIndex::index( entity.type() ), indexSet_.index( entity ) );
  }

public:
  explicit CheckPartitionDataHandle ( const GridView &gridView )
    : gridView_( gridView ),
      rank_( gridView_.comm().rank() ),
      indexSet_( gridView_.indexSet() ),
      invalidDimension_( false ),
      invalidCodimension_( false ),
      invalidEntity_( false ),
      invalidSendEntity_( false ),
      invalidReceiveEntity_( false ),
      invalidSize_( false ),
      selfReceive_( false ),
      doubleInterior_( false ),
      interiorBorder_( false )
  {
    Dune::Hybrid::forEach( std::make_index_sequence< dimension+1 >{},
      [ & ]( auto i ){ contains_[ i ] = Dune::Capabilities::canCommunicate< Grid, i >::v; } );
  }

  ~CheckPartitionDataHandle ()
  {
    bool sendFailure = false;
    bool receiveFailure = false;

    typedef typename GridView::template Codim< 0 >::Iterator Iterator;
    typedef typename GridView::template Codim< 0 >::Entity Entity;
    const Iterator end = gridView_.template end< 0 >();
    for( Iterator it = gridView_.template begin< 0 >(); it != end; ++it )
    {
      const Entity entity = *it;

      if( entity.partitionType() == Dune::InteriorEntity )
        continue;

      const bool wasSent = (sendSet_.find( key( entity ) ) != sendSet_.end());
      if( Dune::EntityCommHelper< iftype >::send( entity.partitionType() ) && !wasSent )
      {
        std::cout << "[ " << rank_ << " ] Error: No data sent on non-interior entity "
                  << grid().globalIdSet().id( entity ) << " of partition type "
                  << entity.partitionType() << " contained in the send set." << std::endl;
        sendFailure = true;
      }

      const bool wasReceived = (receiveSet_.find( key( entity ) ) != receiveSet_.end());
      if( Dune::EntityCommHelper< iftype >::receive( entity.partitionType() ) && !wasReceived )
      {
        std::cout << "[ " << rank_ << " ] Error: No data received on non-interior entity "
                  << grid().globalIdSet().id( entity ) << " of partition type "
                  << entity.partitionType() << " contained in the receive set." << std::endl;
        receiveFailure = true;
      }
    }

    if( invalidDimension_ )
      std::cerr << "[ " << rank_ << " ] Error: Invalid dimension passed during communication." << std::endl;
    if( invalidCodimension_ )
      std::cerr << "[ " << rank_ << " ] Error: Invalid codimension passed during communication." << std::endl;
    if( invalidEntity_ )
      std::cerr << "[ " << rank_ << " ] Error: Uncontained entity passed during communication." << std::endl;
    if( invalidSendEntity_ )
      std::cerr << "[ " << rank_ << " ] Error: Sent data on entity not contained in communication interface." << std::endl;
    if( invalidReceiveEntity_ )
      std::cerr << "[ " << rank_ << " ] Error: Received data on entity not contained in communication interface." << std::endl;
    if( invalidSize_ )
      std::cerr << "[ " << rank_ << " ] Error: Wrong size passed during communication." << std::endl;
    if( selfReceive_ )
      std::cerr << "[ " << rank_ << " ] Warning: Received data from own process during communication." << std::endl;
    if( doubleInterior_ )
      std::cerr << "[ " << rank_ << " ] Error: Received interior data on interior entity." << std::endl;
    if( interiorBorder_ )
      std::cerr << "[ " << rank_ << " ] Error: Received interior data on border entity / border data on interior entity." << std::endl;
    if( sendFailure )
      std::cerr << "[ " << rank_ << " ] Error: No data sent on a non-interior entity within the send set." << std::endl;
    if( receiveFailure )
      std::cerr << "[ " << rank_ << " ] Error: No data received on a non-interior entity within the receive set." << std::endl;
  }

  bool contains ( const int dim, const int codim ) const
  {
    invalidDimension_ |= (dim != dimension);
    invalidCodimension_ |= ((codim < 0) || (codim > dimension));
    return ((codim >= 0) && (codim <= dimension) ? contains_[ codim ] : false);
  }

  bool fixedSize ( const int dim, const int codim ) const
  {
    invalidDimension_ |= (dim != dimension);
    invalidCodimension_ |= ((codim < 0) || (codim > dimension));
    return true;
  }

  template< class Entity >
  std::size_t size ([[maybe_unused]] const Entity &entity) const
  {
    static_assert( (Entity::dimension == dimension), "Entity has invalid dimension." );
    static_assert( (Entity::codimension >= 0) || (Entity::codimension <= dimension), "Entity has invalid codimension." );
    return (contains_[ Entity::codimension ] ? 2 : 0);
  }

  template< class Buffer, class Entity >
  void gather ( Buffer &buffer, const Entity &entity ) const
  {
    static_assert( (Entity::dimension == dimension), "Entity has invalid dimension." );
    static_assert( (Entity::codimension >= 0) || (Entity::codimension <= dimension), "Entity has invalid codimension." );

    invalidEntity_ |= !contains_[ Entity::codimension ];
    invalidSendEntity_ |= !Dune::EntityCommHelper< iftype >::send( entity.partitionType() );

    buffer.write( rank_ );
    buffer.write( int( entity.partitionType() ) );

    sendSet_.insert( key( entity ) );
  }

  template< class Buffer, class Entity >
  void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
  {
    static_assert( (Entity::dimension == dimension), "Entity has invalid dimension." );
    static_assert( (Entity::codimension >= 0) || (Entity::codimension <= dimension), "Entity has invalid codimension." );

    invalidEntity_ |= !contains_[ Entity::codimension ];
    invalidSize_ |= (n != size( entity ));
    invalidReceiveEntity_ |= !Dune::EntityCommHelper< iftype >::receive( entity.partitionType() );

    int rank, partitionType;
    buffer.read( rank );
    buffer.read( partitionType );

    receiveSet_.insert( key( entity ) );

    selfReceive_ |= (rank == rank_);
    if( (partitionType == int( Dune::InteriorEntity )) && (entity.partitionType() == Dune::InteriorEntity) )
    {
      std::cout << "[ " << rank_ << " ] Error: Receive interior data from process " << rank
                << " on interior entity " << grid().globalIdSet().id( entity ) << "." << std::endl;
      doubleInterior_ = true;
    }
    if( (partitionType == int( Dune::BorderEntity )) && (entity.partitionType() == Dune::InteriorEntity) )
    {
      std::cout << "[ " << rank_ << " ] Error: Receive border data from process " << rank
                << " on interior entity " << grid().globalIdSet().id( entity ) << "." << std::endl;
      interiorBorder_ = true;
    }
    if( (partitionType == int( Dune::InteriorEntity )) && (entity.partitionType() == Dune::BorderEntity) )
    {
      std::cout << "[ " << rank_ << " ] Error: Receive interior data from process " << rank
                << " on border entity " << grid().globalIdSet().id( entity ) << "." << std::endl;
      interiorBorder_ = true;
    }
  }

  static void apply ( const GridView &gridView )
  {
    std::cout << "Checking communication for " << iftype << "..." << std::endl;
    CheckPartitionDataHandle handle( gridView );
    auto commFuture = gridView.communicate( handle, iftype, Dune::ForwardCommunication );
    if( ! commFuture.ready() )
      commFuture.wait();
  }

  const Grid &grid () const { return gridView_.grid(); }

private:
  const GridView &gridView_;
  const int rank_;
  const IndexSet &indexSet_;
  std::bitset< dimension+1 > contains_;
  mutable CommSet sendSet_, receiveSet_;
  mutable bool invalidDimension_;
  mutable bool invalidCodimension_;
  mutable bool invalidEntity_;
  mutable bool invalidSendEntity_;
  bool invalidReceiveEntity_;
  bool invalidSize_;
  bool selfReceive_;
  bool doubleInterior_;
  bool interiorBorder_;
};



/** \brief check the partition types of the periodic copies of an entity on this process
 *
 * An interior entity is owned by this process and does not lie on a process
 * boundary. Hence, all its periodic copies on this process are overlap, front,
 * or ghost entities. A border entity can have border copies, e.g., the
 * vertices on opposite sides of a periodic domain.
 *
 * Only local coincidences are checked; copies on other processes are not visible here.
 */
template< class GridView >
inline void checkPeriodicPartitionTypes ( const GridView &gridView,
                                          const Dune::GridCheck::PeriodicMap< GridView > &periodicMap )
{
  Dune::Hybrid::forEach( std::make_index_sequence< GridView::dimension+1 >{}, [ & ]( auto codim )
  {
    if constexpr( Dune::Capabilities::hasEntity< typename GridView::Grid, codim >::v )
    {
      typedef Dune::GridCheck::PeriodicIdentification< GridView > Identification;
      const Identification periodic( gridView, codim, periodicMap );
      if( !periodic.hasCoincidences() )
        return;

      std::map< typename Identification::Key, Dune::PartitionType > partitionType;
      for( const auto &element : elements( gridView ) )
        for( unsigned int i = 0; i < element.subEntities( codim ); ++i )
          partitionType[ periodic.subKey( element, i, codim ) ] = element.template subEntity< codim >( i ).partitionType();

      for( std::size_t c = 0; c < periodic.size(); ++c )
      {
        int interior = 0, border = 0;
        for( const auto &key : periodic.members( c ) )
        {
          interior += (partitionType[ key ] == Dune::InteriorEntity);
          border += (partitionType[ key ] == Dune::BorderEntity);
        }
        if( (interior > 1) || ((interior == 1) && (border > 0)) )
          DUNE_THROW( Dune::GridError, "Codim " << codim << " entity has " << interior << " interior and "
                      << border << " border periodic copies (expected: at most one interior and no border copy"
                      << " or no interior copy)." );
      }
    }
  } );
}


/** \brief check the partition types and the communication interfaces of a grid view
 *
 * \param periodicMap maps each coordinate to its representative in the periodic
 *                    unit cell, see dune/grid/test/periodicmap.hh
 */
template< class GridView >
inline void checkPartitionType ( const GridView &gridView,
                                 const Dune::GridCheck::PeriodicMap< GridView > &periodicMap = Dune::GridCheck::noPeriodicMap< GridView >() )
{
  checkPeriodicPartitionTypes( gridView, periodicMap );

  CheckPartitionType< GridView, Dune::Interior_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::InteriorBorder_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::Overlap_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::OverlapFront_Partition >::apply( gridView );
  CheckPartitionType< GridView, Dune::Ghost_Partition >::apply( gridView );

  CheckPartitionDataHandle< GridView, Dune::InteriorBorder_InteriorBorder_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::InteriorBorder_All_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::Overlap_OverlapFront_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::Overlap_All_Interface >::apply( gridView );
  CheckPartitionDataHandle< GridView, Dune::All_All_Interface >::apply( gridView );
}

#endif // DUNE_GRID_TEST_CHECKPARTITION_HH
