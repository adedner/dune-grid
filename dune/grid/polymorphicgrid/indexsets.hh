// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRID_INDEXSETS_HH
#define DUNE_VIRTUALIZEDGRID_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the Polymorphic class
 */

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/polymorphicgrid/idtype.hh>
#include <dune/grid/polymorphicgrid/common/typeerasure.hh>

#include <vector>

namespace Dune {

  template<class GridImp>
  class PolymorphicIndexSet;

  template<class GridImp>
  struct PolymorphicIndexSetDefinition
  {
    using Types = typename IndexSet<GridImp, PolymorphicIndexSet<GridImp>>::Types;
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
      virtual InterfaceImpl *clone () const = 0;
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

    template<class Derived, class I, int codim>
    struct ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using WrappedIndexSet = std::decay_t<I>;

      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using WrappedEntity = typename WrappedIndexSet::template Codim<codim>::Entity;

      int index (Codim<codim>, const Entity& e) const final
      {
        return derived().impl().index(Polymorphic::asWrapped<WrappedEntity>(e));
      }

      int subIndex (Codim<codim>, const Entity& e, int i, int cd) const final
      {
        return derived().impl().template subIndex<codim>(Polymorphic::asWrapped<WrappedEntity>(e), i, cd);
      }

      bool contains (Codim<codim>, const Entity& e) const final
      {
        return derived().impl().contains(Polymorphic::asWrapped<WrappedEntity>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct ImplementationImpl
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      ImplementationImpl (I&& i)
        : impl_( std::forward<I>(i) )
      {}

      ImplementationImpl *clone() const override { return new ImplementationImpl(*this); }

      std::size_t size (int codim) const final
      {
        return impl().size(codim);
      }

      std::size_t size (GeometryType type) const final
      {
        return impl().size(type);
      }

      Types types (int codim) const final
      {
        auto t = impl().types(codim);
        return Types(std::begin(t), std::end(t));
      }

      const auto& impl () const { return impl_; }
      auto& impl () { return impl_; }

    private:
      I impl_;
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
  };

  /** \todo Take the index types from the host grid */
  template<class GridImp>
  class PolymorphicIndexSet :
    public IndexSet<GridImp, PolymorphicIndexSet<GridImp>>
  {
    using Definition = PolymorphicIndexSetDefinition<GridImp>;
    template <class I>
    using Implementation = typename Definition::template Implementation<I>;
    using Interface = typename Definition::Interface;

  public:
    using Types = typename IndexSet<GridImp, PolymorphicIndexSet<GridImp>>::Types;

    enum { dim = GridImp::dimension };

  public:
    template<class Impl>
    explicit PolymorphicIndexSet (Impl&& impl)
      : impl_(new Implementation<Impl>(std::forward<Impl>(impl)))
    {}

    PolymorphicIndexSet (const PolymorphicIndexSet& other)
      : impl_(other.impl_ ? other.impl_->clone() : nullptr)
    {}

    PolymorphicIndexSet (PolymorphicIndexSet &&) = default;

    PolymorphicIndexSet& operator= (const PolymorphicIndexSet& other)
    {
      impl_.reset(other.impl_ ? other.impl_->clone() : nullptr);
      return *this;
    }


    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      return impl().index(Codim<codim>{}, e);
    }

    //! get index of subEntity of a codim 0 entity
    template<int codim>
    int subIndex (const typename GridImp::Traits::template Codim<codim>::Entity& e, int i, int cd) const
    {
      return impl().subIndex(Codim<codim>{}, e, i, cd);
    }

    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const
    {
      return impl().size(codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (GeometryType type) const
    {
      return impl().size(type);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      return impl().types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      static constexpr int codim = EntityType::codimension;
      return impl().contains(Codim<codim>{}, e);
    }

    Interface& impl() const
    {
      return *impl_;
    }

  private:
    std::unique_ptr<Interface> impl_;
  };


  template <class GridImp>
  struct PolymorphicIdSetDefinition
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
      virtual InterfaceImpl *clone () const = 0;
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


    template<class Derived, class I, int codim>
    struct ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using WrappedIdSet = std::decay_t<I>;
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using WrappedEntity = typename WrappedIdSet::template Codim<codim>::Entity;

      IdType id (Codim<codim>, const Entity& e) const final
      {
        return derived().impl().template id<codim>(Polymorphic::asWrapped<WrappedEntity>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      using WrappedIdSet = std::decay_t<I>;
      using Entity = typename GridImp::Traits::template Codim<0>::Entity;
      using WrappedEntity = typename WrappedIdSet::template Codim<0>::Entity;

      ImplementationImpl (I&& i)
        : impl_( std::forward<I>(i) )
      {}

      ImplementationImpl* clone () const final { return new ImplementationImpl(*this); }

      IdType subId (const Entity& e, int i, int codim) const final
      {
        return impl().subId(Polymorphic::asWrapped<WrappedEntity>(e), i, codim);
      }

      const auto& impl () const { return impl_; }
      auto& impl () { return impl_; }

    private:
      I impl_;
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
  };

  template <class GridImp>
  class PolymorphicIdSet :
    public IdSet<GridImp,PolymorphicIdSet<GridImp>,PolymorphicIdType>
  {
    using Definition = PolymorphicIdSetDefinition<GridImp>;
    template <class I>
    using Implementation = typename Definition::template Implementation<I>;
    using Interface = typename Definition::Interface;

  public:
    using IdType = typename GridImp::Traits::GlobalIdSet::IdType;

    enum { dim = GridImp::dimension };

  public:
    template<class Impl>
    explicit PolymorphicIdSet(Impl&& impl)
      : impl_( new Implementation<Impl>( std::forward<Impl>(impl) ) )
    {}

    PolymorphicIdSet(const PolymorphicIdSet& other)
      : impl_(other.impl_ ? other.impl_->clone() : nullptr)
    {}

    PolymorphicIdSet (PolymorphicIdSet&&) = default;

    PolymorphicIdSet& operator= (const PolymorphicIdSet& other)
    {
      impl_.reset(other.impl_ ? other.impl_->clone() : nullptr);
      return *this;
    }


    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      return impl().id(Codim<cd>{}, e);
    }

    //! get id of subEntity
    IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      return impl().subId(e, i, codim);
    }

    Interface& impl() const
    {
      return *impl_;
    }

  private:
    std::unique_ptr<Interface> impl_;
  };

}  // namespace Dune


#endif
