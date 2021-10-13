// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_ENTITYOWNER_HH
#define DUNE_GRID_UTILITY_ENTITYOWNER_HH

#include <algorithm>
#include <cassert>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>

namespace Dune {

/** \brief Attributes used in the generic overlap model
 *
   \code
    #include <dune/grid/common/gridenums.hh>
    \endcode
  *
  * The values are ordered intentionally in order to be able to
  * define ranges of partition types.
  *
  * @ingroup GIRelatedTypes
  */
enum class OwnerType {
  Owner=0,          //!< all interior and owned border entities
  Overlap=1,        //!< all entities lying in the overlap zone and non-owned border entities
  Ghost=2,          //!< all ghost and front entities
  Unknown=10        //!< No unique owner type can be specified
};


/// Calculate unique owner rank for all entities in a given GridView
template<class GridView>
class EntityOwner
{
private:
  /// A DataHandle class to calculate the minimum of a vector which is accompanied by an index set
  /**
   * \tparam IS  An indexSet of the GridView
   * \tparam V   The vector type to compute the elementwise minimum
   **/
  template<class IS, class Vec>
  class MinimumExchange
      : public Dune::CommDataHandleIF<MinimumExchange<IS,Vec>,typename Vec::value_type>
  {
  public:
    /** \brief export type of data for message buffer */
    using ValueType = typename Vec::value_type;

    /** \brief constructor. */
    MinimumExchange (const IS& indexset, Vec& vec)
      : indexset_(indexset)
      , vec_(vec)
    {}

    /** \brief returns true if data for this codim should be communicated */
    bool contains (int /*dim*/, int /*codim*/) const
    {
      return true;
    }

    /** \brief returns true if size per entity of given dim and codim is a constant */
    bool fixedSize (int /*dim*/, int /*codim*/) const
    {
      return true ;
    }

    /** \brief number of values to send */
    template <class Entity>
    std::size_t size (const Entity& e) const
    {
      return 1 ;
    }

    /** \brief pack data from user to message buffer */
    template <class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& e) const
    {
      buff.write(vec_[indexset_.index(e)]);
    }

    /** \brief Unpack data from message buffer to user
     *
     * \param n The number of objects sent by the sender
     */
    template <class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& e, [[maybe_unused]] std::size_t n)
    {
      assert(n == 1);
      ValueType x;
      buff.read(x);
      if (x >= 0) // other is -1 means, he does not want it
        vec_[indexset_.index(e)] = std::min(x, vec_[indexset_.index(e)]);
    }

  private:
    const IS& indexset_;
    Vec& vec_;
  };

  static MCMGLayout layout ()
  {
    return [](GeometryType /*gt*/, int /*dim*/) { return true; };
  }

public:

  using IndexSet = MultipleCodimMultipleGeomTypeMapper<GridView>;

public:
  /** \brief Constructor constructs a mapper for all entities */
  EntityOwner (const GridView& gridView)
    : indexSet_{gridView, layout()}
    , assignment_(indexSet_.size(), -1)
    , rank_{gridView.comm().rank()}
  {
    // assign own rank to entities that I might have
    for (auto const& e : elements(gridView))
    {
      Dune::Hybrid::forEach(Dune::StaticIntegralRange<int,GridView::dimension+1>{}, [&](auto cd)
      {
        for (unsigned int i = 0; i < e.subEntities(cd); ++i)
        {
          Dune::PartitionType subPartitionType = e.template subEntity<cd>(i).partitionType();

          assignment_[indexSet_.subIndex(e,i,cd)]
            = (subPartitionType == Dune::PartitionType::InteriorEntity ||
               subPartitionType == Dune::PartitionType::BorderEntity)
            ? rank_   // set to own rank
            : - 1;    // it is a ghost entity, I will not possibly own it.
        }
      });
    }

    // exchange entity rank through communication
    MinimumExchange dh{indexSet_, assignment_};
    gridView.communicate(dh,
      Dune::InterfaceType::InteriorBorder_All_Interface,
      Dune::CommunicationDirection::ForwardCommunication);
  }

  /** \brief Which rank is the entity assigned to? */
  template <class Entity>
  int owner (const Entity& entity) const
  {
    return assignment_[indexSet_.index(entity)];
  }

  /** \brief Which rank is the `i`-th subentity of `entity` of given `codim` assigned to? */
  template <class Entity>
  int owner (const Entity& entity, int i, unsigned int codim) const
  {
    return assignment_[indexSet_.subIndex(entity,i,codim)];
  }


  /** \brief OwnerType of the given `entity`. */
  template <class Entity>
  Dune::OwnerType ownerType (const Entity& entity) const
  {
    Dune::PartitionType pt = entity.partitionType();
    if (pt != Dune::PartitionType::BorderEntity)
      return partitionTypeToOwnerType(pt);
    else
      return owner(entity) == rank_
        ? Dune::OwnerType::Owner
        : Dune::OwnerType::Overlap;
  }

  /** \brief OwnerType of the `i`-th subentity of `entity` of given `codim`. */
  template <class Entity>
  Dune::OwnerType ownerType (const Entity& entity, int i, unsigned int codim) const
  {
    static_assert(Entity::codimension == 0);
    Dune::PartitionType pt = Dune::Hybrid::switchCases(
      std::make_index_sequence<GridView::dimension+1>{}, codim,
      [&](auto cd) { return entity.template subEntity<cd>(i).partitionType(); },
      [&] {
          DUNE_THROW(Dune::Exception, "Invalid codimension codim=" << codim);
          return Dune::PartitionType{};
      });

    if (pt != Dune::PartitionType::BorderEntity)
      return partitionTypeToOwnerType(pt);
    else
      return owner(entity,i,codim) == rank_
        ? Dune::OwnerType::Owner
        : Dune::OwnerType::Overlap;
  }

private:
  static Dune::OwnerType partitionTypeToOwnerType (Dune::PartitionType pt)
  {
    switch (pt) {
      case Dune::PartitionType::InteriorEntity:
        return Dune::OwnerType::Owner;
      case Dune::PartitionType::OverlapEntity:
        return Dune::OwnerType::Overlap;
      case Dune::PartitionType::FrontEntity:
      case Dune::PartitionType::GhostEntity:
        return Dune::OwnerType::Ghost;
      default:
        DUNE_THROW(Dune::Exception,
          "PartitionType " << pt << " cannot uniquely be mapped to OwnerType.");
        return Dune::OwnerType::Unknown;
    }
  }

private:
  IndexSet indexSet_;
  std::vector<int> assignment_;
  int rank_;
};

} // end namespace Dune

#endif // DUNE_GRID_UTILITY_ENTITYOWNER_HH