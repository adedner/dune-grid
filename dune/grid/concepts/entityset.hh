#ifndef DUNE_GRID_CONCEPTS_ENTITY_SET_HH
#define DUNE_GRID_CONCEPTS_ENTITY_SET_HH

#include <dune/grid/concepts/entity.hh>
#include <dune/grid/concepts/entityiterator.hh>
#include <dune/grid/concepts/intersection.hh>
#include <dune/grid/concepts/intersectioniterator.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/concepts/indexset.hh>

#include <dune/grid/common/gridenums.hh>

namespace Dune::Concept {

namespace Impl {

  template<class ES>
  concept EntitySetEssentials = requires(ES es, Dune::GeometryType type, int codim)
  {
    typename ES::Traits;
    typename ES::ctype;
    requires IndexSet<typename ES::IndexSet>;
    requires Intersection<typename ES::Intersection>;
    requires IntersectionIterator<typename ES::IntersectionIterator>;
    { ES::conforming        } -> std::convertible_to< bool                                  >;
    { ES::dimension         } -> std::convertible_to< int                                   >;
    { ES::dimensionworld    } -> std::convertible_to< int                                   >;
    { es.grid()             } -> std::convertible_to< const typename ES::Grid&              >;
    { es.indexSet()         } -> std::convertible_to< const typename ES::IndexSet&          >;
    { es.size(codim)        } -> std::convertible_to< int                                   >;
    { es.size(type)         } -> std::convertible_to< int                                   >;
    { es.comm()             } -> std::convertible_to< typename ES::CollectiveCommunication  >;
    { es.overlapSize(codim) } -> std::convertible_to< int                                   >;
    { es.ghostSize(codim)   } -> std::convertible_to< int                                   >;
    //! gv.communicate(std::declval<Handler&>(),std::declval<Dune::InterfaceType>(), std::declval<CommunicationDirection>()), // FIXME use a default handler to instantiate this function
    requires std::copy_constructible<ES>;
  };

}

/**
 * @brief Model of an entity set
 * @ingroup GridConcepts
 * @details An entity set is an iterable set of entities whose index set can
 *    contiguously index all the iterated entities. A Dune::GridView is a
 *    template for this model, but not necessarily the only one. In particular,
 *    entity sets are used to enumerate sub-ranges of entities contained in the
 *    grid view.
 */
template<class ES, int codim>
concept EntitySet = requires(ES es, Dune::GeometryType type, const typename ES::template Codim<codim>::Entity& entity)
{
  requires Impl::EntitySetEssentials<ES>;
  requires (codim != 0) || requires {
    { es.ibegin(entity)     } -> std::convertible_to< typename ES::IntersectionIterator>;
    { es.iend(entity)       } -> std::convertible_to< typename ES::IntersectionIterator>;
  };
  requires Geometry<typename ES::template Codim<codim>::Geometry>;
  requires Geometry<typename ES::template Codim<codim>::LocalGeometry>;
  requires EntityIterator<typename ES::template Codim<codim>::Iterator>;
  { es.template begin<codim>()            } -> std::convertible_to< typename ES::template Codim<codim>::Iterator>;
  { es.template end<codim>()              } -> std::convertible_to< typename ES::template Codim<codim>::Iterator>;
};

}  // end namespace Dune::Concept

#endif // DUNE_GRID_CONCEPTS_ENTITY_SET_HH
