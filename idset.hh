// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_GRID_IDSET_HH
#define DUNE_PYTHON_GRID_IDSET_HH

#include <functional>
#include <type_traits>

#include <dune/common/typeutilities.hh>

#include <dune/python/common/string.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    namespace detail
    {

      // GridId
      // ------

      template< class Id >
      struct GridId
      {
        Id id;
      };



      // Comparison Operators
      // --------------------

      template< class Id >
      inline static bool operator== ( const GridId< Id > &a, const GridId< Id > &b )
      {
        return (a.id == b.id);
      }

      template< class Id >
      inline static bool operator!= ( const GridId< Id > &a, const GridId< Id > &b )
      {
        return !(a == b);
      }

      template< class Id >
      inline static bool operator< ( const GridId< Id > &a, const GridId< Id > &b )
      {
        return (a.id < b.id);
      }

      template< class Id >
      inline static bool operator>= ( const GridId< Id > &a, const GridId< Id > &b )
      {
        return !(a.id < b.id);
      }

      template< class Id >
      inline static auto operator<= ( const GridId< Id > &a, const GridId< Id > &b )
        -> std::enable_if_t< std::is_same< decltype( a.id <= b.id ), bool >::value, bool >
      {
        return (a.id <= b.id);
      }

      template< class Id >
      inline static auto operator<= ( const GridId< Id > &a, const GridId< Id > &b )
        -> std::enable_if_t< !std::is_same< decltype( a.id <= b.id ), bool >::value, bool >
      {
        return (a.id < b.id) || (a.id == b.id);
      }

      template< class Id >
      inline static bool operator> ( const GridId< Id > &a, const GridId< Id > &b )
      {
        return !(a.id <= b.id);
      }



      // to_string
      // ---------

      template< class Id >
      inline static auto to_string ( const GridId< Id > &id, PriorityTag< 2 > )
        -> decltype( static_cast< std::string >( id.id ) )
      {
        return static_cast< std::string >( id.id );
      }

      template< class Id >
      inline static auto to_string ( const GridId< Id > &id, PriorityTag< 1 > )
        -> decltype( to_string( id.id ) )
      {
        return to_string( id.id );
      }

      template< class Id >
      inline static std::string to_string ( const GridId< Id > &id, PriorityTag< 0 > )
      {
        std::ostringstream s;
        s << id.id;
        return s.str();
      }

      template< class Id >
      inline static std::string to_string ( const GridId< Id > &id )
      {
        return to_string( id, PriorityTag< 42 >() );
      }



      // registerSubId
      // -------------

      template< class Entity, class IdSet, class... options >
      inline static auto registerSubId ( pybind11::class_< IdSet, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< Entity::codimension == 0 >
      {
        cls.def( "subId", [] ( const IdSet &idSet, const Entity &entity, int i, int codim ) {
            if( (codim < Entity::codimension) || (codim > Entity::dimension) )
              throw pybind11::value_error( "Invalid codimension: " + std::to_string( codim ) + " (must be in [" + std::to_string( Entity::codimension ) + ", " + std::to_string( Entity::dimension ) + "])" );
            return idSet.subId( entity, i, codim );
          } );
      }

      template< class Entity, class IdSet, class... options >
      inline static void registerSubId ( pybind11::class_< IdSet, options... >, PriorityTag< 0 > )
      {}

      template< class Entity, class IdSet, class... options >
      inline static void registerSubId ( pybind11::class_< IdSet, options... > cls )
      {
        return registerSubId< Entity >( cls, PriorityTag< 42 >() );
      }

    } // namespace detail



    // registerGridIdSet
    // -----------------

    template< class Grid, class IdSet, class... options >
    inline static void registerGridIdSet ( pybind11::handle scope, pybind11::class_< IdSet, options... > cls )
    {
      typedef detail::GridId< typename IdSet::IdType > Id;
      auto id = insertClass< Id >( cls, "Id", GenerateTypeName( "Dune::Python::detail::GridId", GenerateTypeName( cls, "IdType" ) ) );
      if( id.second )
      {
        id.first.def( pybind11::self == pybind11::self );
        id.first.def( pybind11::self != pybind11::self );

        id.first.def( pybind11::self < pybind11::self );
        id.first.def( pybind11::self <= pybind11::self );
        id.first.def( pybind11::self > pybind11::self );
        id.first.def( pybind11::self >= pybind11::self );

        std::hash< typename IdSet::IdType > hash;
        id.first.def( "__hash__", [ hash ] ( const Id &self ) { return pybind11::int_( hash( self.id ) ); } );
        id.first.def( "__str__", [] ( const Id &self ) { return to_string( self ); } );
      }

      Hybrid::forEach( std::make_integer_sequence< int, Grid::dimension+1 >(), [ &cls ] ( auto &&codim ) {
          typedef typename Grid::template Codim< codim >::Entity Entity;

          using pybind11::operator""_a;

          cls.def( "id", [] ( const IdSet &self, const Entity &entity ) { return self.id( entity ); }, "entity"_a );
          detail::registerSubId< Entity >( cls );
        } );
    }



    // registerGridIdSets
    // ------------------

    template< class Grid, class... options >
    inline static void registerHierarchicalGridIdSets ( pybind11::class_< Grid, options... > cls )
    {
      typedef typename Grid::LocalIdSet LocalIdSet;
      auto local = insertClass< LocalIdSet >( cls, "LocalIdSet", GenerateTypeName( cls, "LocalIdSet" ) );
      if( local.second )
        registerGridIdSet< Grid >( cls, local.first );
      cls.def_property_readonly( "localIdSet", [] ( const Grid &self ) -> const LocalIdSet & { return self.localIdSet(); }, pybind11::keep_alive< 0, 1 >() );

      typedef typename Grid::GlobalIdSet GlobalIdSet;
      auto global = insertClass< GlobalIdSet >( cls, "GlobalIdSet", GenerateTypeName( cls, "GlobalIdSet" ) );
      if( global.second )
        registerGridIdSet< Grid >( cls, global.first );
      cls.def_property_readonly( "globalIdSet", [] ( const Grid &self ) -> const GlobalIdSet & { return self.globalIdSet(); }, pybind11::keep_alive< 0, 1 >() );
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_IDSET_HH
