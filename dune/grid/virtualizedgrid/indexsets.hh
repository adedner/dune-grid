// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRID_INDEXSETS_HH
#define DUNE_VIRTUALIZEDGRID_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the VirtualizedGrid class
 */

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/virtualizedgrid/idtype.hh>
#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

#include <vector>

namespace Dune {

  template<class GridImp>
  class VirtualizedGridIndexSet;

  template<class GridImp>
  struct VirtualizedGridIndexSetDefinition
  {
    using Types = typename IndexSet<GridImp, VirtualizedGridIndexSet<GridImp>>::Types;
    enum {dim = GridImp::dimension};

    template<int codim>
    struct InterfaceCodim
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;

      virtual ~InterfaceCodim () = default;
      virtual int index (Codim<codim>, const Entity& e) const = 0;
      virtual int subIndex (Codim<codim>, const Entity& e, int i, int cd) const = 0;
      virtual bool contains (Codim<codim>, const Entity& e) const = 0;
    };

    template<int... codims>
    struct InterfaceImpl
        : virtual InterfaceCodim<codims>...
    {
      virtual ~InterfaceImpl () = default;
      virtual std::size_t size (int codim) const = 0;
      virtual std::size_t size (GeometryType type) const = 0;
      virtual Types types (int codim) const = 0;

      using InterfaceCodim<codims>::index...;
      using InterfaceCodim<codims>::subIndex...;
      using InterfaceCodim<codims>::contains...;
    };

    template<class Seq>
    struct Interface_t;

    template<int... codims>
    struct Interface_t<std::integer_sequence<int,codims...>>
    {
      using type = InterfaceImpl<codims...>;
    };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;

    template<class Derived, class Wrapper, int codim>
    struct ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using WrappedIndexSet = typename Wrapper::Wrapped;

      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using WrappedEntity = typename WrappedIndexSet::template Codim<codim>::Entity;

      int index (Codim<codim>, const Entity& e) const final
      {
        return derived().get().index(Polymorphic::asWrapped<WrappedEntity>(e));
      }

      int subIndex (Codim<codim>, const Entity& e, int i, int cd) const final
      {
        return derived().get().template subIndex<codim>(Polymorphic::asWrapped<WrappedEntity>(e), i, cd);
      }

      bool contains (Codim<codim>, const Entity& e) const final
      {
        return derived().get().contains(Polymorphic::asWrapped<WrappedEntity>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Wrapper, int... codims>
    struct ImplementationImpl
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<Wrapper,codims...>, Wrapper, codims>...
      , public Wrapper
    {
      using Wrapper::Wrapper;

      std::size_t size (int codim) const final
      {
        return this->get().size(codim);
      }

      std::size_t size (GeometryType type) const final
      {
        return this->get().size(type);
      }

      Types types (int codim) const final
      {
        auto t = this->get().types(codim);
        return Types(std::begin(t), std::end(t));
      }
    };

    template<class I, class Seq>
    struct Implementation_t;

    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>>
    {
      using type = ImplementationImpl<I,codims...>;
    };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dim+1>>::type;

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };

  /** \todo Take the index types from the host grid */
  template<class GridImp>
  class VirtualizedGridIndexSet :
    public VirtualizedGridIndexSetDefinition<GridImp>::Base,
    public IndexSet<GridImp, VirtualizedGridIndexSet<GridImp>>
  {
    using Definition = VirtualizedGridIndexSetDefinition<GridImp>;
    using Base = typename Definition::Base;

  public:
    using Types = typename IndexSet<GridImp, VirtualizedGridIndexSet<GridImp>>::Types;

    enum { dim = GridImp::dimension };

  public:
    template <class Impl, disableCopyMove<VirtualizedGridIndexSet,Impl> = 0>
    VirtualizedGridIndexSet (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      return this->asInterface().index(Codim<codim>{}, e);
    }

    //! get index of subEntity of a codim 0 entity
    template<int codim>
    int subIndex (const typename GridImp::Traits::template Codim<codim>::Entity& e, int i, int cd) const
    {
      return this->asInterface().subIndex(Codim<codim>{}, e, i, cd);
    }

    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const
    {
      return this->asInterface().size(codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (GeometryType type) const
    {
      return this->asInterface().size(type);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      return this->asInterface().types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      static constexpr int codim = EntityType::codimension;
      return this->asInterface().contains(Codim<codim>{}, e);
    }
  };


  template <class GridImp>
  struct VirtualizedGridIdSetDefinition
  {
    using IdType = typename GridImp::Traits::GlobalIdSet::IdType;

    enum { dim = GridImp::dimension };

    template<int codim>
    struct InterfaceCodim
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;

      virtual ~InterfaceCodim () = default;
      virtual IdType id (Codim<codim>, const Entity& e) const = 0;
    };

    template<int... codims>
    struct InterfaceImpl
        : virtual InterfaceCodim<codims>...
    {
      using Entity = typename GridImp::Traits::template Codim<0>::Entity;

      virtual ~InterfaceImpl () = default;
      virtual IdType subId (const Entity& e, int i, int codim) const = 0;

      using InterfaceCodim<codims>::id...;
    };

    template<class Seq>
    struct Interface_t;

    template<int... codims>
    struct Interface_t<std::integer_sequence<int,codims...>>
    {
      using type = InterfaceImpl<codims...>;
    };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;


    template<class Derived, class Wrapper, int codim>
    struct ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using WrappedIdSet = typename Wrapper::Wrapped;

      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using WrappedEntity = typename WrappedIdSet::template Codim<codim>::Entity;

      IdType id (Codim<codim>, const Entity& e) const final
      {
        return derived().get().template id<codim>(Polymorphic::asWrapped<WrappedEntity>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class Wrapper, int... codims>
    struct ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<Wrapper,codims...>, Wrapper, codims>...
      , public Wrapper
    {
      using WrappedIdSet = typename Wrapper::Wrapped;

      using Entity = typename GridImp::Traits::template Codim<0>::Entity;
      using WrappedEntity = typename WrappedIdSet::template Codim<0>::Entity;

      using Wrapper::Wrapper;

      IdType subId (const Entity& e, int i, int codim) const final
      {
        return this->get().subId(Polymorphic::asWrapped<WrappedEntity>(e), i, codim);
      }
    };

    template<class I, class Seq>
    struct Implementation_t;

    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>>
    {
      using type = ImplementationImpl<I,codims...>;
    };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dim+1>>::type;

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };

  template <class GridImp>
  class VirtualizedGridIdSet :
    public VirtualizedGridIdSetDefinition<GridImp>::Base,
    public IdSet<GridImp,VirtualizedGridIdSet<GridImp>,VirtualizedGridIdType>
  {
    using Definition = VirtualizedGridIdSetDefinition<GridImp>;
    using Base = typename Definition::Base;

  public:
    using IdType = typename GridImp::Traits::GlobalIdSet::IdType;

    enum { dim = GridImp::dimension };

  public:
    template <class Impl, disableCopyMove<VirtualizedGridIdSet,Impl> = 0>
    VirtualizedGridIdSet (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      return this->asInterface().id(Codim<cd>{}, e);
    }

    //! get id of subEntity
    IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      return this->asInterface().subId(e, i, codim);
    }
  };

}  // namespace Dune


#endif
