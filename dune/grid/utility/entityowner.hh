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
#include <dune/grid/utility/ownertype.hh>

namespace Dune {

/// Calculate unique owner rank for all entities in a given GridView
template<class GridView>
class EntityOwner
{
private:
  using IndexSet = MultipleCodimMultipleGeomTypeMapper<GridView>;
  using Container = std::vector<int>;

  /// A DataHandle class to calculate the minimum of a vector which is accompanied by an index set
  class MinimumExchange
      : public Dune::CommDataHandleIF<MinimumExchange,int>
  {
  public:
    /** \brief constructor. */
    MinimumExchange (const IndexSet& indexSet, Container& ranks, int codim, int commSize)
      : indexSet_{indexSet}
      , ranks_{ranks}
      , codim_{codim}
      , commSize_{commSize}
    {}

    /** \brief returns true if data for this codim should be communicated */
    bool contains (int /*dim*/, int codim) const
    {
      return codim_ < 0 || codim_ == codim;
    }

    /** \brief returns true if size per entity of given dim and codim is a constant */
    bool fixedSize (int /*dim*/, int /*codim*/) const
    {
      return true;
    }

    /** \brief number of values to send */
    template <class Entity>
    std::size_t size (const Entity& e) const
    {
      return 1;
    }

    /** \brief pack data from user to message buffer */
    template <class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& e) const
    {
      assert(codim_ < 0 || Entity::codimension == codim_);
      buff.write(ranks_[indexSet_.index(e)]);
    }

    /** \brief Unpack data from message buffer to user
     *
     * \param n The number of objects sent by the sender
     */
    template <class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& e, [[maybe_unused]] std::size_t n)
    {
      assert(n == 1);
      assert(codim_ < 0 || Entity::codimension == codim_);

      int otherRank = -1;
      buff.read(otherRank);
      assert(otherRank >= 0);

      auto const idx = indexSet_.index(e);
      ranks_[idx] = std::min(otherRank, ranks_[idx]);
    }

  private:
    const IndexSet& indexSet_;
    Container& ranks_;
    int codim_;
    int commSize_;
  };

  static MCMGLayout makeLayout (int codim)
  {
    return [codim](GeometryType gt, int dim) { return codim < 0 || int(dim - gt.dim()) == codim; };
  }


public:
  /** \brief Constructor constructs a mapper for all entities (if codim < 0) or entities of the
   *  specified codimension `codim` only.
   */
  EntityOwner (const GridView& gridView, int codim = -1)
    : indexSet_{gridView, makeLayout(codim)}
    , assignment_(indexSet_.size(), gridView.comm().size())
    , codim_{codim}
    , rank_{gridView.comm().rank()}
    , size_{gridView.comm().size()}
  {
    auto const layout = makeLayout(codim_);
    for (auto const& e : elements(gridView))
    {
      auto refElem = referenceElement(e);
      Hybrid::forEach(StaticIntegralRange<int,GridView::dimension+1>{}, [&](auto cd) {
        for (int i = 0; i < refElem.size(cd); ++i) {
          if (layout(refElem.type(i,cd), GridView::dimension)) {
            PartitionType subPartitionType = e.template subEntity<cd>(i).partitionType();

            assignment_[indexSet_.subIndex(e,i,cd)]
              = (subPartitionType == PartitionType::InteriorEntity ||
                 subPartitionType == PartitionType::BorderEntity)
              ? rank_       // set to own rank
              : size_;      // it is a ghost/overlap entity, thus I will not own it.
          }
        }
      });
    }

    // exchange entity rank through communication
    MinimumExchange dh{indexSet_, assignment_, codim_, size_};
    gridView.communicate(dh,
      Dune::InterfaceType::InteriorBorder_All_Interface,
      Dune::CommunicationDirection::ForwardCommunication);
  }


  /** \brief Which rank is the entity assigned to? */
  template <class Entity>
  int ownerRank (const Entity& entity) const
  {
    assert(codim_ < 0 || Entity::codimension == codim_);
    return assignment_[indexSet_.index(entity)];
  }

  /** \brief Which rank is the `i`-th subentity of `entity` of given `codim` assigned to? */
  template <class Entity>
  int ownerRank (const Entity& entity, int i, unsigned int codim) const
  {
    assert(codim_ < 0 || int(Entity::codimension + codim) == codim_);
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
      return ownerRank(entity) == rank_
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
      return ownerRank(entity,i,codim) == rank_
        ? Dune::OwnerType::Owner
        : Dune::OwnerType::Overlap;
  }

  /** \brief Return number of entities (in the layout) owned by the current rank */
  std::size_t countOwner () const
  {
    return std::count(assignment_.begin(), assignment_.end(), rank_);
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
  Container assignment_;
  int codim_;
  int rank_;
  int size_;
};

} // end namespace Dune

#endif // DUNE_GRID_UTILITY_ENTITYOWNER_HH