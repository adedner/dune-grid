// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_VIRTUALIZEDGRID_GRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_GRID_HH

/** \file
 * \brief The VirtualizedGrid class
 */

#include <string>
#include <map>
#include <any>

#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/grid.hh>

// The components of the VirtualizedGrid interface
#include "datahandle.hh"
#include "entity.hh"
#include "entityiterator.hh"
#include "entityseed.hh"
#include "geometry.hh"
#include "gridview.hh"
#include "idtype.hh"
#include "indexsets.hh"
#include "intersectioniterator.hh"

#if HAVE_MPI
  #include <dune/common/parallel/mpicommunication.hh>
  using VirtualizedCommunication = Dune::Communication<MPI_Comm>;
#else
  using VirtualizedCommunication = Dune::Communication<No_Comm>;
#endif

namespace Dune
{
  // Forward declaration
  template<int dimension, int dimensionworld, typename ct = double>
  class VirtualizedGrid;

  template<int dimension, int dimensionworld, typename ct>
  struct VirtualizedGridFamily;

  template<class Grid, int dimension, int dimensionworld, typename ct>
  class VirtualizedGridDefinition
  {
    template <PartitionIteratorType p>
    struct _Partition {};

    template <PartitionIteratorType... pp>
    struct _Partitions {};

    template <class... TT>
    struct _Types {};

  public:
    //! type of the used GridFamily for this grid
    using GridFamily = VirtualizedGridFamily<dimension, dimensionworld, ct>;

    //! the Traits
    using Traits = typename GridFamily::Traits;

    template<class DataType>
    struct InterfaceDataType
    {
      virtual ~InterfaceDataType () = default;
      virtual void communicate (VirtualizedCommDataHandle<DataType,Grid>&, InterfaceType iftype, CommunicationDirection dir) const = 0;
      virtual void communicate (VirtualizedCommDataHandle<DataType,Grid>&, InterfaceType iftype, CommunicationDirection dir, int level) const = 0;
    };

    template<class DataTypes>
    struct InterfaceDataTypes;

    template<class... DataTypes>
    struct InterfaceDataTypes<_Types<DataTypes...>>
        : virtual InterfaceDataType<DataTypes>...
    {
      virtual ~InterfaceDataTypes () = default;
      using InterfaceDataType<DataTypes>::communicate...;
    };

    template<int codim, PartitionIteratorType pitype>
    struct InterfaceCodimPartition
    {
      virtual ~InterfaceCodimPartition () = default;

      using LevelIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator;
      using LeafIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator;

      virtual LevelIterator lbegin (Codim<codim>, _Partition<pitype>, int level) const = 0;
      virtual LevelIterator lend (Codim<codim>, _Partition<pitype>, int level) const = 0;
      virtual LeafIterator leafbegin (Codim<codim>, _Partition<pitype>) const = 0;
      virtual LeafIterator leafend (Codim<codim>, _Partition<pitype>) const = 0;
    };

    template<int codim, class Partitions>
    struct InterfaceCodim;

    template<int codim, PartitionIteratorType... pitypes>
    struct InterfaceCodim<codim, _Partitions<pitypes...>>
        : virtual InterfaceCodimPartition<codim,pitypes>...
    {
      virtual ~InterfaceCodim () = default;

      using LevelIterator = typename Traits::template Codim<codim>::LevelIterator;
      using LeafIterator = typename Traits::template Codim<codim>::LeafIterator;
      using Entity = typename Traits::template Codim<codim>::Entity;
      using EntitySeed = typename Traits::template Codim<codim>::EntitySeed;

      virtual LevelIterator lbegin (Codim<codim>, int level) const = 0;
      virtual LevelIterator lend (Codim<codim>, int level) const = 0;
      virtual LeafIterator leafbegin (Codim<codim>) const = 0;
      virtual LeafIterator leafend (Codim<codim>) const = 0;
      virtual Entity entity (Codim<codim>, const EntitySeed& seed) const = 0;

      using InterfaceCodimPartition<codim,pitypes>::lbegin...;
      using InterfaceCodimPartition<codim,pitypes>::lend...;
      using InterfaceCodimPartition<codim,pitypes>::leafbegin...;
      using InterfaceCodimPartition<codim,pitypes>::leafend...;
    };

    using AllPartitions = _Partitions<Interior_Partition, InteriorBorder_Partition, Overlap_Partition, OverlapFront_Partition, All_Partition, Ghost_Partition>;

    using AllDataTypes = _Types<std::byte,char,unsigned char,signed char,short,unsigned short,int,unsigned int,long,unsigned long,long long,unsigned long long,float,double,long double>;

    template<int... codims>
    struct InterfaceImpl
        : virtual InterfaceCodim<codims,AllPartitions>...
        , virtual InterfaceDataTypes<AllDataTypes>
    {
      using Entity0 = typename Traits::template Codim<0>::Entity;

      virtual ~InterfaceImpl () = default;
      virtual int maxLevel () const = 0;
      virtual int size (int level, int codim) const = 0;
      virtual int size (int codim) const = 0;
      virtual int size (int level, GeometryType type) const = 0;
      virtual int size (GeometryType type) const = 0;
      virtual size_t numBoundarySegments () const = 0;
      virtual const typename Traits::GlobalIdSet& globalIdSet () const = 0;
      virtual const typename Traits::LocalIdSet& localIdSet () const = 0;
      virtual const typename Traits::LevelIndexSet& levelIndexSet (int level) const = 0;
      virtual const typename Traits::LeafIndexSet& leafIndexSet () const = 0;
      virtual void globalRefine (int refCount) = 0;
      virtual bool mark (int refCount, const typename Traits::template Codim<0>::Entity & e) = 0;
      virtual int getMark (const typename Traits::template Codim<0>::Entity & e) const = 0;
      virtual bool preAdapt () = 0;
      virtual bool adapt () = 0;
      virtual void postAdapt () = 0;
      virtual unsigned int overlapSize (int codim) const = 0;
      virtual unsigned int ghostSize (int codim) const = 0;
      virtual unsigned int overlapSize (int level, int codim) const = 0;
      virtual unsigned int ghostSize (int level, int codim) const = 0;
      virtual const VirtualizedCommunication& comm () const = 0;

      virtual typename Traits::LevelIntersectionIterator ilevelbegin (const Entity0& entity) const = 0;
      virtual typename Traits::LevelIntersectionIterator ilevelend (const Entity0& entity) const = 0;
      virtual typename Traits::LeafIntersectionIterator ileafbegin (const Entity0& entity) const = 0;
      virtual typename Traits::LeafIntersectionIterator ileafend (const Entity0& entity) const = 0;

      using InterfaceCodim<codims,AllPartitions>::lbegin...;
      using InterfaceCodim<codims,AllPartitions>::lend...;
      using InterfaceCodim<codims,AllPartitions>::leafbegin...;
      using InterfaceCodim<codims,AllPartitions>::leafend...;
      using InterfaceCodim<codims,AllPartitions>::entity...;
      using InterfaceDataTypes<AllDataTypes>::communicate;
    };

    template<class Derived, class DataType>
    struct ImplementationDataType
        : virtual InterfaceDataType<DataType>
    {
      void communicate (VirtualizedCommDataHandle<DataType,Grid>& dh, InterfaceType iftype, CommunicationDirection dir) const final
      {
        derived().get().leafGridView().communicate(dh,iftype,dir);
      }

      void communicate (VirtualizedCommDataHandle<DataType,Grid>& dh, InterfaceType iftype, CommunicationDirection dir, int level) const final
      {
        derived().get().levelGridView(level).communicate(dh,iftype,dir);
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Derived, class DataTypes>
    struct ImplementationDataTypes;

    template<class Derived, class... DataTypes>
    struct ImplementationDataTypes<Derived,_Types<DataTypes...>>
        : public ImplementationDataType<Derived,DataTypes>...
    {};

    template<class Derived, class WrappedGrid, int codim, PartitionIteratorType pitype>
    struct ImplementationCodimPartition
        : virtual InterfaceCodimPartition<codim,pitype>
    {
      using LevelIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator;
      using LevelIteratorImpl = typename LevelIterator::Implementation;
      using LeafIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator;
      using LeafIteratorImpl = typename LeafIterator::Implementation;

      virtual LevelIterator lbegin (Codim<codim>, _Partition<pitype>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImpl{derived().get().levelGridView(level).template begin<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      virtual LevelIterator lend (Codim<codim>, _Partition<pitype>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImpl{derived().get().levelGridView(level).template end<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      virtual LeafIterator leafbegin (Codim<codim>, _Partition<pitype>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImpl{derived().get().leafGridView().template begin<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

      virtual LeafIterator leafend (Codim<codim>, _Partition<pitype>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImpl{derived().get().leafGridView().template end<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Derived, class WrappedGrid, int codim, class Partitions>
    struct ImplementationCodim;

    template<class Derived, class WrappedGrid, int codim, PartitionIteratorType... pitypes>
    struct ImplementationCodim<Derived,WrappedGrid,codim,_Partitions<pitypes...>>
        : virtual InterfaceCodim<codim,_Partitions<pitypes...>>
        , public ImplementationCodimPartition<Derived,WrappedGrid,codim,pitypes>...
    {
      using LevelIterator = typename Traits::template Codim<codim>::LevelIterator;
      using LevelIteratorImpl = typename LevelIterator::Implementation;
      using LeafIterator = typename Traits::template Codim<codim>::LeafIterator;
      using LeafIteratorImpl = typename LeafIterator::Implementation;
      using Entity = typename Traits::template Codim<codim>::Entity;
      using EntityImpl = typename Entity::Implementation;
      using EntitySeed = typename Traits::template Codim<codim>::EntitySeed;
      using EntitySeedImpl = typename EntitySeed::Implementation;

      virtual LevelIterator lbegin (Codim<codim>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImpl{derived().get().levelGridView(level).template begin<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      virtual LevelIterator lend (Codim<codim>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImpl{derived().get().levelGridView(level).template end<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      virtual LeafIterator leafbegin (Codim<codim>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImpl{derived().get().leafGridView().template begin<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

      virtual LeafIterator leafend (Codim<codim>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImpl{derived().get().leafGridView().template end<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

      virtual Entity entity (Codim<codim>, const EntitySeed& seed) const final
      {
        return EntityImpl{derived().get().entity(Polymorphic::asWrapped<EntitySeedImpl>(seed))};
      }

      using ImplementationCodimPartition<Derived,WrappedGrid,codim,pitypes>::lbegin...;
      using ImplementationCodimPartition<Derived,WrappedGrid,codim,pitypes>::lend...;
      using ImplementationCodimPartition<Derived,WrappedGrid,codim,pitypes>::leafbegin...;
      using ImplementationCodimPartition<Derived,WrappedGrid,codim,pitypes>::leafend...;

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Wrapper, int... codims>
    struct ImplementationImpl
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<typename Wrapper::Wrapped,codims...>, typename Wrapper::Wrapped, codims, AllPartitions>...
      , public ImplementationDataTypes<ImplementationImpl<typename Wrapper::Wrapped,codims...>, AllDataTypes>
      , public Wrapper
    {
      using WrappedGrid = typename Wrapper::Wrapped;
      using Entity0 = typename Traits::template Codim<0>::Entity;
      using WrappedEntity0 = typename WrappedGrid::template Codim<0>::Entity;

      using LevelIntersectionIterator = typename Traits::LevelIntersectionIterator;
      using LevelIntersectionIteratorImpl = typename LevelIntersectionIterator::Implementation;
      using LeafIntersectionIterator = typename Traits::LeafIntersectionIterator;
      using LeafIntersectionIteratorImpl = typename LeafIntersectionIterator::Implementation;

      template<class TT, disableCopyMove<ImplementationImpl, TT> = 0>
      ImplementationImpl (TT&& tt)
        : Wrapper{std::forward<TT>(tt)}
        , globalIdSet_( Wrapper::get().globalIdSet() )
        , localIdSet_( Wrapper::get().localIdSet() )
        , leafIndexSet_( Wrapper::get().leafIndexSet() )
        , comm_( Wrapper::get().comm() )
      {
        for (int i = 0; i <= maxLevel(); i++)
        {
          levelIndexSets_.push_back(
            new VirtualizedGridIndexSet<const Grid>( this->get().levelIndexSet(i) )
          );
        }
      }

      ~ImplementationImpl () noexcept
      {
        for (size_t i = 0; i < levelIndexSets_.size(); i++)
          if (levelIndexSets_[i])
            delete (levelIndexSets_[i]);
      }

      virtual int maxLevel () const final { return this->get().maxLevel(); }

      virtual LevelIntersectionIterator ilevelbegin (const Entity0& entity) const final
      {
        return LevelIntersectionIteratorImpl{this->get().levelGridView(entity.level()).ibegin(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      virtual LevelIntersectionIterator ilevelend (const Entity0& entity) const final
      {
        return LevelIntersectionIteratorImpl{this->get().levelGridView(entity.level()).iend(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      virtual LeafIntersectionIterator ileafbegin (const Entity0& entity) const final
      {
        return LeafIntersectionIteratorImpl{this->get().leafGridView().ibegin(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      virtual LeafIntersectionIterator ileafend (const Entity0& entity) const final
      {
        return LeafIntersectionIteratorImpl{this->get().leafGridView().iend(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      virtual int size (int level, int codim) const final { return this->get().size(level, codim); }
      virtual int size (int codim) const final { return this->get().size(codim); }
      virtual int size (int level, GeometryType type) const final { return this->get().size(level, type); }
      virtual int size (GeometryType type) const final { return this->get().size(type); }
      virtual size_t numBoundarySegments () const final { return this->get().numBoundarySegments(); }

      virtual const typename Traits::GlobalIdSet& globalIdSet () const final
      {
        return dynamic_cast<const typename Traits::GlobalIdSet&>(globalIdSet_);
      }

      virtual const typename Traits::LocalIdSet& localIdSet () const final
      {
        return dynamic_cast<const typename Traits::LocalIdSet&>(localIdSet_);
      }

      virtual const typename Traits::LevelIndexSet& levelIndexSet (int level) const final
      {
        return dynamic_cast<const typename Traits::LevelIndexSet&>(*levelIndexSets_[level]);
      }

      virtual const typename Traits::LeafIndexSet& leafIndexSet () const final
      {
        return dynamic_cast<const typename Traits::LeafIndexSet&>(leafIndexSet_);
      }

      virtual void globalRefine (int refCount) final { return this->get().globalRefine(refCount); }

      virtual bool mark (int refCount, const Entity0& e) final
      {
        return this->get().mark(refCount, Polymorphic::asWrapped<WrappedEntity0>(e));
      }

      virtual int getMark (const Entity0 & e) const final
      {
        return this->get().getMark(Polymorphic::asWrapped<WrappedEntity0>(e));
      }

      virtual bool preAdapt () final { return this->get().preAdapt(); }
      virtual bool adapt () final { return this->get().adapt(); }
      virtual void postAdapt () final { return this->get().postAdapt(); }
      virtual unsigned int overlapSize (int codim) const final { return this->get().leafGridView().overlapSize(codim); }
      virtual unsigned int ghostSize (int codim) const final { return this->get().leafGridView().ghostSize(codim); }
      virtual unsigned int overlapSize (int level, int codim) const final { return this->get().levelGridView(level).overlapSize(codim); }
      virtual unsigned int ghostSize (int level, int codim) const final { return this->get().levelGridView(level).ghostSize(codim); }
      virtual const VirtualizedCommunication& comm () const override { return comm_; }

    private:
      VirtualizedGridIdSet<const Grid> globalIdSet_;
      VirtualizedGridIdSet<const Grid> localIdSet_;
      std::vector<VirtualizedGridIndexSet<const Grid>*> levelIndexSets_;
      VirtualizedGridIndexSet<const Grid> leafIndexSet_;
      VirtualizedCommunication comm_{};
    };


    template<class Seq>
    struct Interface_t;

    template<int... codims>
    struct Interface_t<std::integer_sequence<int,codims...>>
    {
      using type = InterfaceImpl<codims...>;
    };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dimension+1>>::type;

    template<class I, class Seq>
    struct Implementation_t;

    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>>
    {
      using type = ImplementationImpl<I,codims...>;
    };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dimension+1>>::type;

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };

} // namespace Dune

#endif // DUNE_GRID_VIRTUALIZEDGRID_GRID_HH
