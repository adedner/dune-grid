// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_GRID_VIRTUALIZEDGRID_GRID_HH
#define DUNE_GRID_VIRTUALIZEDGRID_GRID_HH

/** \file
 * \brief The Polymorphic class
 */

#include <string>
#include <map>
#include <any>

#include <dune/common/parallel/communication.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/grid.hh>

// The components of the Polymorphic interface
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
  using PolymorphicCommunication = Dune::Communication<MPI_Comm>;
#else
  using PolymorphicCommunication = Dune::Communication<No_Comm>;
#endif

namespace Dune
{
  // Forward declaration
  template<int dimension, int dimensionworld, typename ct = double>
  class PolymorphicGrid;

  template<int dimension, int dimensionworld, typename ct>
  struct PolymorphicGridFamily;

  template<class Grid, int dimension, int dimensionworld, typename ct>
  struct PolymorphicGridDefinition
  {
    template <PartitionIteratorType p>
    struct _Partition {};

    template <PartitionIteratorType... pp>
    struct _Partitions {};

    template <class... TT>
    struct _Types {};

  public:
    //! type of the used GridFamily for this grid
    using GridFamily = PolymorphicGridFamily<dimension, dimensionworld, ct>;

    //! the Traits
    using Traits = typename GridFamily::Traits;

    template<class DataType>
    struct InterfaceDataType
    {
      virtual ~InterfaceDataType () = default;
      virtual void communicate (PolymorphicCommDataHandle<DataType,Grid>&, InterfaceType iftype, CommunicationDirection dir) const = 0;
      virtual void communicate (PolymorphicCommDataHandle<DataType,Grid>&, InterfaceType iftype, CommunicationDirection dir, int level) const = 0;
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
      virtual InterfaceImpl *clone () const = 0;
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
      virtual const PolymorphicCommunication& comm () const = 0;

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

    template<class Derived, class I, class DataType>
    struct ImplementationDataType
        : virtual InterfaceDataType<DataType>
    {
      using DataHandleImp = PolymorphicCommDataHandle<DataType,Grid>;

      void communicate (DataHandleImp& dh, InterfaceType iftype, CommunicationDirection dir) const final
      {
        derived().impl().leafGridView().communicate(dh,iftype,dir);
      }

      void communicate (DataHandleImp& dh, InterfaceType iftype, CommunicationDirection dir, int level) const final
      {
        derived().impl().levelGridView(level).communicate(dh,iftype,dir);
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Derived, class I, class DataTypes>
    struct ImplementationDataTypes;

    template<class Derived, class I, class... DataTypes>
    struct ImplementationDataTypes<Derived,I,_Types<DataTypes...>>
        : public ImplementationDataType<Derived,I,DataTypes>...
    {};

    template<class Derived, class I, int codim, PartitionIteratorType pitype>
    struct ImplementationCodimPartition
        : virtual InterfaceCodimPartition<codim,pitype>
    {
      using WrappedGrid = std::decay_t<I>;
      using LevelIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LevelIterator;
      using LevelIteratorImp = typename LevelIterator::Implementation;
      using LeafIterator = typename Traits::template Codim<codim>::template Partition<pitype>::LeafIterator;
      using LeafIteratorImp = typename LeafIterator::Implementation;

      LevelIterator lbegin (Codim<codim>, _Partition<pitype>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImp{derived().impl().levelGridView(level).template begin<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      LevelIterator lend (Codim<codim>, _Partition<pitype>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImp{derived().impl().levelGridView(level).template end<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      LeafIterator leafbegin (Codim<codim>, _Partition<pitype>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImp{derived().impl().leafGridView().template begin<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

      LeafIterator leafend (Codim<codim>, _Partition<pitype>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImp{derived().impl().leafGridView().template end<codim,pitype>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Derived, class I, int codim, class Partitions>
    struct ImplementationCodim;

    template<class Derived, class I, int codim, PartitionIteratorType... pitypes>
    struct ImplementationCodim<Derived,I,codim,_Partitions<pitypes...>>
        : virtual InterfaceCodim<codim,_Partitions<pitypes...>>
        , public ImplementationCodimPartition<Derived,I,codim,pitypes>...
    {
      using WrappedGrid = std::decay_t<I>;
      using LevelIterator = typename Traits::template Codim<codim>::LevelIterator;
      using LevelIteratorImp = typename LevelIterator::Implementation;
      using LeafIterator = typename Traits::template Codim<codim>::LeafIterator;
      using LeafIteratorImp = typename LeafIterator::Implementation;
      using Entity = typename Traits::template Codim<codim>::Entity;
      using EntityImpl = typename Entity::Implementation;
      using EntitySeed = typename Traits::template Codim<codim>::EntitySeed;
      using WrappedEntitySeed = typename WrappedGrid::Traits::template Codim<codim>::EntitySeed;

      LevelIterator lbegin (Codim<codim>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImp{derived().impl().levelGridView(level).template begin<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      LevelIterator lend (Codim<codim>, int level) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LevelIteratorImp{derived().impl().levelGridView(level).template end<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LevelIterator{};
        }
      }

      LeafIterator leafbegin (Codim<codim>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImp{derived().impl().leafGridView().template begin<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

      LeafIterator leafend (Codim<codim>) const final
      {
        if constexpr(Dune::Capabilities::hasEntityIterator<WrappedGrid,codim>::v)
          return LeafIteratorImp{derived().impl().leafGridView().template end<codim>()};
        else {
          DUNE_THROW(Dune::NotImplemented, "EntityIterator<codim="<<codim<<"> not implemented");
          return LeafIterator{};
        }
      }

      Entity entity (Codim<codim>, const EntitySeed& seed) const final
      {
        return EntityImpl{derived().impl().entity(Polymorphic::asWrapped<WrappedEntitySeed>(seed))};
      }

      using ImplementationCodimPartition<Derived,I,codim,pitypes>::lbegin...;
      using ImplementationCodimPartition<Derived,I,codim,pitypes>::lend...;
      using ImplementationCodimPartition<Derived,I,codim,pitypes>::leafbegin...;
      using ImplementationCodimPartition<Derived,I,codim,pitypes>::leafend...;

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct ImplementationImpl
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims, AllPartitions>...
      , public ImplementationDataTypes<ImplementationImpl<I,codims...>, I, AllDataTypes>
    {
      using WrappedGrid = std::decay_t<I>;
      using Entity0 = typename Traits::template Codim<0>::Entity;
      using WrappedEntity0 = typename WrappedGrid::template Codim<0>::Entity;

      using LevelIntersectionIterator = typename Traits::LevelIntersectionIterator;
      using LevelIntersectionIteratorImp = typename LevelIntersectionIterator::Implementation;
      using LeafIntersectionIterator = typename Traits::LeafIntersectionIterator;
      using LeafIntersectionIteratorImp = typename LeafIntersectionIterator::Implementation;

      ImplementationImpl (I&& i)
        : impl_(std::forward<I>(i))
        , globalIdSet_( impl().globalIdSet() )
        , localIdSet_( impl().localIdSet() )
        , leafIndexSet_( impl().leafIndexSet() )
        , comm_( impl().comm() )
      {
        for (int i = 0; i <= maxLevel(); i++)
        {
          PolymorphicIndexSet<const Grid>* p
            = new PolymorphicIndexSet<const Grid>( impl().levelIndexSet(i) );
          levelIndexSets_.push_back(p);
        }
      }

      ~ImplementationImpl ()
      {
        for (size_t i = 0; i < levelIndexSets_.size(); i++)
          if (levelIndexSets_[i])
            delete (levelIndexSets_[i]);
      }

      ImplementationImpl* clone () const final { return new ImplementationImpl( *this ); }
      int maxLevel () const final { return impl().maxLevel(); }

      LevelIntersectionIterator ilevelbegin (const Entity0& entity) const final
      {
        return LevelIntersectionIteratorImp{impl()
          .levelGridView(entity.level())
          .ibegin(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      LevelIntersectionIterator ilevelend (const Entity0& entity) const final
      {
        return LevelIntersectionIteratorImp{impl()
          .levelGridView(entity.level())
          .iend(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      LeafIntersectionIterator ileafbegin (const Entity0& entity) const final
      {
        return LeafIntersectionIteratorImp{impl()
          .leafGridView()
          .ibegin(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      LeafIntersectionIterator ileafend (const Entity0& entity) const final
      {
        return LeafIntersectionIteratorImp{impl()
          .leafGridView()
          .iend(Polymorphic::asWrapped<WrappedEntity0>(entity))};
      }

      int size (int level, int codim) const final { return impl().size(level, codim); }
      int size (int codim) const final { return impl().size(codim); }
      int size (int level, GeometryType type) const final { return impl().size(level, type); }
      int size (GeometryType type) const final { return impl().size(type); }
      size_t numBoundarySegments () const final { return impl().numBoundarySegments(); }

      const typename Traits::GlobalIdSet& globalIdSet () const final
      {
        return dynamic_cast<const typename Traits::GlobalIdSet&>(globalIdSet_);
      }

      const typename Traits::LocalIdSet& localIdSet () const final
      {
        return dynamic_cast<const typename Traits::LocalIdSet&>(localIdSet_);
      }

      const typename Traits::LevelIndexSet& levelIndexSet (int level) const final
      {
        return dynamic_cast<const typename Traits::LevelIndexSet&>(*levelIndexSets_[level]);
      }

      const typename Traits::LeafIndexSet& leafIndexSet () const final
      {
        return dynamic_cast<const typename Traits::LeafIndexSet&>(leafIndexSet_);
      }

      void globalRefine (int refCount) final { return impl().globalRefine(refCount); }

      bool mark (int refCount, const Entity0& entity) final
      {
        return impl().mark(refCount, Polymorphic::asWrapped<WrappedEntity0>(entity));
      }

      int getMark (const Entity0& entity) const final
      {
        return impl().getMark(Polymorphic::asWrapped<WrappedEntity0>(entity));
      }

      bool preAdapt () final { return impl().preAdapt(); }
      bool adapt () final { return impl().adapt(); }
      void postAdapt () final { return impl().postAdapt(); }
      unsigned int overlapSize (int codim) const final { return impl().leafGridView().overlapSize(codim); }
      unsigned int ghostSize (int codim) const final { return impl().leafGridView().ghostSize(codim); }
      unsigned int overlapSize (int level, int codim) const final { return impl().levelGridView(level).overlapSize(codim); }
      unsigned int ghostSize (int level, int codim) const final { return impl().levelGridView(level).ghostSize(codim); }
      const PolymorphicCommunication& comm () const override { return comm_; }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
      PolymorphicIdSet<const Grid> globalIdSet_;
      PolymorphicIdSet<const Grid> localIdSet_;
      std::vector<PolymorphicIndexSet<const Grid>*> levelIndexSets_;
      PolymorphicIndexSet<const Grid> leafIndexSet_;
      PolymorphicCommunication comm_{};
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
  };

} // namespace Dune

#endif // DUNE_GRID_VIRTUALIZEDGRID_GRID_HH
