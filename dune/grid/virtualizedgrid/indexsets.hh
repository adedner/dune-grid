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

#include <vector>

namespace Dune {

  /** \todo Take the index types from the host grid */
  template<class GridImp>
  class VirtualizedGridLevelIndexSet :
    public IndexSet<GridImp, VirtualizedGridLevelIndexSet<GridImp>>
  {
  public:

    typedef typename IndexSet<GridImp, VirtualizedGridLevelIndexSet<GridImp>>::Types Types;

    enum {dim = GridImp::dimension};

  private:
    // VIRTUALIZATION BEGIN
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
    struct Interface_t<std::integer_sequence<int,codims...>> { using type = InterfaceImpl<codims...>; };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;

    template<class Derived, class I, int codim>
    struct DUNE_PRIVATE ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<codim,dim,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<codim>::Entity>;

      int index (Codim<codim>, const Entity& e) const final {
        return derived().impl().index(upcast<EntityImpl>(e));
      }
      int subIndex (Codim<codim>, const Entity& e, int i, int cd) const final {
        return derived().impl().template subIndex<codim>(upcast<EntityImpl>(e), i, cd);
      }
      bool contains (Codim<codim>, const Entity& e) const final {
        return derived().impl().contains(upcast<EntityImpl>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct DUNE_PRIVATE ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      ImplementationImpl ( I&& i ) : impl_( std::forward<I>(i) ) {}
      ImplementationImpl *clone() const override { return new ImplementationImpl( *this ); }

      std::size_t size (int codim) const override { return impl().size(codim); }
      std::size_t size (GeometryType type) const override { return impl().size(type); }
      Types types (int codim) const override { return impl().types(codim); }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };

    template<class I, class Seq>
    struct Implementation_t;
    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>> { using type = ImplementationImpl<I,codims...>; };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dim+1>>::type;
    // VIRTUALIZATION END

  public:
    template< class ImplLevelIndexSet >
    explicit VirtualizedGridLevelIndexSet(ImplLevelIndexSet&& implLevelIndexSet)
    : impl_( new Implementation<ImplLevelIndexSet>( std::forward<ImplLevelIndexSet>( implLevelIndexSet ) ) )
    {}

    VirtualizedGridLevelIndexSet(const VirtualizedGridLevelIndexSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLevelIndexSet ( VirtualizedGridLevelIndexSet && ) = default;

    VirtualizedGridLevelIndexSet& operator=(const VirtualizedGridLevelIndexSet& other) {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const {
      return impl_->index(Codim<codim>{}, e);
    }

    //! get index of subEntity of a codim 0 entity
    template<int codim>
    int subIndex (const typename GridImp::Traits::template Codim<codim>::Entity& e, int i, int cd) const {
      return impl_->subIndex(Codim<codim>{}, e, i, cd);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const {
      return impl_->size(codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (GeometryType type) const {
      return impl_->size(type);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const {
      return impl_->types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const {
      static constexpr int codim = EntityType::codimension;
      return impl_->contains(Codim<codim>{}, e);
    }

    std::unique_ptr<Interface> impl_;
  };


  template<class GridImp>
  class VirtualizedGridLeafIndexSet :
    public IndexSet<GridImp, VirtualizedGridLeafIndexSet<GridImp>>
  {

  public:
    typedef typename IndexSet<GridImp, VirtualizedGridLeafIndexSet<GridImp>>::Types Types;

    enum {dim = GridImp::dimension};

  private:
    // VIRTUALIZATION BEGIN
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
    struct Interface_t<std::integer_sequence<int,codims...>> { using type = InterfaceImpl<codims...>; };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;

    template<class Derived, class I, int codim>
    struct DUNE_PRIVATE ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<codim,dim,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<codim>::Entity>;

      int index (Codim<codim>, const Entity& e) const final {
        return derived().impl().index(upcast<EntityImpl>(e));
      }
      int subIndex (Codim<codim>, const Entity& e, int i, int cd) const final {
        return derived().impl().template subIndex<codim>(upcast<EntityImpl>(e), i, cd);
      }
      bool contains (Codim<codim>, const Entity& e) const final {
        return derived().impl().contains(upcast<EntityImpl>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct DUNE_PRIVATE ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      ImplementationImpl ( I&& i ) : impl_( std::forward<I>(i) ) {}
      ImplementationImpl *clone() const override { return new ImplementationImpl( *this ); }

      std::size_t size (int codim) const override { return impl().size(codim); }
      std::size_t size (GeometryType type) const override { return impl().size(type); }
      Types types (int codim) const override { return impl().types(codim); }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };

    template<class I, class Seq>
    struct Implementation_t;
    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>> { using type = ImplementationImpl<I,codims...>; };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dim+1>>::type;
    // VIRTUALIZATION END


  public:
    template< class ImplLeafIndexSet >
    explicit VirtualizedGridLeafIndexSet(ImplLeafIndexSet&& implLeafIndexSet)
    : impl_( new Implementation<ImplLeafIndexSet>( std::forward<ImplLeafIndexSet>(implLeafIndexSet) ) )
    {}

    VirtualizedGridLeafIndexSet(const VirtualizedGridLeafIndexSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLeafIndexSet ( VirtualizedGridLeafIndexSet && ) = default;

    VirtualizedGridLeafIndexSet& operator=(const VirtualizedGridLeafIndexSet& other) {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }


    //! get index of an entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int codim>
    int index (const typename GridImp::template Codim<codim>::Entity& e) const {
      return impl_->index(Codim<codim>{}, e);
    }


    //! get index of subEntity of a codim 0 entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int codim>
    int subIndex (const typename GridImp::Traits::template Codim<codim>::Entity& e, int i, int cd) const {
      return impl_->subIndex(Codim<codim>{}, e, i, cd);
    }


    //! get number of entities of given type
    std::size_t size (GeometryType type) const {
      return impl_->size(type);
    }


    //! get number of entities of given codim
    std::size_t size (int codim) const {
      return impl_->size(codim);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const {
      return impl_->types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      static constexpr int codim = EntityType::codimension;
      return impl_->contains(Codim<codim>{}, e);
    }

    std::unique_ptr<Interface> impl_;
  };




  template <class GridImp>
  class VirtualizedGridGlobalIdSet :
    public IdSet<GridImp,VirtualizedGridGlobalIdSet<GridImp>,
        typename GridImp::Traits::LocalIdSet::IdType>
  {
  public:
    //! define the type used for persistent indices
    typedef typename GridImp::Traits::GlobalIdSet::IdType IdType;

    enum {dim = GridImp::dimension};

  private:
    // VIRTUALIZATION BEGIN
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
    struct Interface_t<std::integer_sequence<int,codims...>> { using type = InterfaceImpl<codims...>; };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;


    template<class Derived, class I, int codim>
    struct DUNE_PRIVATE ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<codim,dim,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<codim>::Entity>;

      IdType id (Codim<codim>, const Entity& e) const final {
        return derived().impl().template id<codim>(upcast<EntityImpl>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct DUNE_PRIVATE ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      using Entity = typename GridImp::Traits::template Codim<0>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<0,dim,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<0>::Entity>;

      ImplementationImpl ( I&& i ) : impl_( std::forward<I>(i) ) {}
      ImplementationImpl *clone() const override { return new ImplementationImpl( *this ); }

      IdType subId (const Entity& e, int i, int codim) const override {
        return impl().subId(upcast<EntityImpl>(e), i, codim);
      }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };

    template<class I, class Seq>
    struct Implementation_t;
    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>> { using type = ImplementationImpl<I,codims...>; };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dim+1>>::type;
    // VIRTUALIZATION END


  public:
    template< class ImplGlobalIdSet >
    explicit VirtualizedGridGlobalIdSet(ImplGlobalIdSet&& implGlobalIdSet)
    : impl_( new Implementation<ImplGlobalIdSet>( std::forward<ImplGlobalIdSet>( implGlobalIdSet ) ) )
    {}

    VirtualizedGridGlobalIdSet(const VirtualizedGridGlobalIdSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridGlobalIdSet ( VirtualizedGridGlobalIdSet && ) = default;

    VirtualizedGridGlobalIdSet& operator= (const VirtualizedGridGlobalIdSet& other) {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const {
      return impl_->id(Codim<cd>{}, e);
    }

    //! get id of subEntity
    IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const {
      return impl_->subId(e, i, codim);
    }

    std::unique_ptr<Interface> impl_;
  };




  template<class GridImp>
  class VirtualizedGridLocalIdSet :
    public IdSet<GridImp,VirtualizedGridLocalIdSet<GridImp>,
        typename GridImp::Traits::LocalIdSet::IdType>
  {
  public:
    //! define the type used for persistent local ids
    typedef typename GridImp::Traits::LocalIdSet::IdType IdType;

    enum {dim = GridImp::dimension};

  private:
    // VIRTUALIZATION BEGIN
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
    struct Interface_t<std::integer_sequence<int,codims...>> { using type = InterfaceImpl<codims...>; };

    using Interface = typename Interface_t<std::make_integer_sequence<int,dim+1>>::type;


    template<class Derived, class I, int codim>
    struct DUNE_PRIVATE ImplementationCodim
      : virtual InterfaceCodim<codim>
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<codim,dim,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<codim>::Entity>;

      IdType id (Codim<codim>, const Entity& e) const final {
        return derived().impl().template id<codim>(upcast<EntityImpl>(e));
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };

    template<class I, int... codims>
    struct DUNE_PRIVATE ImplementationImpl final
      : virtual InterfaceImpl<codims...>
      , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      using Entity = typename GridImp::Traits::template Codim<0>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<0,dim,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<0>::Entity>;

      ImplementationImpl ( I&& i ) : impl_( std::forward<I>(i) ) {}
      ImplementationImpl *clone() const override { return new ImplementationImpl( *this ); }

      IdType subId (const Entity& e, int i, int codim) const override {
        return impl().subId(upcast<EntityImpl>(e), i, codim);
      }

      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

    private:
      I impl_;
    };

    template<class I, class Seq>
    struct Implementation_t;
    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>> { using type = ImplementationImpl<I,codims...>; };

    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,dim+1>>::type;
    // VIRTUALIZATION END


  public:
    template< class ImplLocalIdSet >
    explicit VirtualizedGridLocalIdSet(ImplLocalIdSet&& implLocalIdSet)
    : impl_( new Implementation<ImplLocalIdSet>( std::forward<ImplLocalIdSet>(implLocalIdSet) ) )
    {}

    VirtualizedGridLocalIdSet(const VirtualizedGridLocalIdSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLocalIdSet ( VirtualizedGridLocalIdSet && ) = default;

    VirtualizedGridLocalIdSet& operator= (const VirtualizedGridLocalIdSet& other) {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const {
      return impl_->id(Codim<cd>{}, e);
    }

    //! get id of subEntity
    IdType subId (const typename GridImp::template Codim<0>::Entity& e, int i, int codim) const {
      return impl_->subId(e, i, codim);
    }

    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune


#endif
