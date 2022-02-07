// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual int index0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual int index1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const = 0;
      virtual int indexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const = 0;
      virtual int subIndex0 (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
      virtual int subIndex1 (const typename GridImp::Traits::template Codim<1>::Entity& e, int i, int codim) const = 0;
      virtual int subIndexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e, int i, int codim) const = 0;
      virtual std::size_t size (int codim) const = 0;
      virtual std::size_t size (GeometryType type) const = 0;
      virtual Types types (int codim) const = 0;
      virtual bool contains0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual bool contains1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const = 0;
      virtual bool containsDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename VirtualizedGridEntity<0, dim, GridImp>::template Implementation<typename I::template Codim<0>::Entity> ImplEntity0;
      typedef typename VirtualizedGridEntity<1, dim, GridImp>::template Implementation<typename I::template Codim<1>::Entity> ImplEntity1;
      typedef typename VirtualizedGridEntity<dim, dim, GridImp>::template Implementation<typename I::template Codim<dim>::Entity> ImplEntityDim;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual int index0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().index(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl()
        );
      }

      virtual int index1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const override
      {
        return impl().index(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl()
        );
      }

      virtual int indexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const override
      {
        return impl().index(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl()
        );
      }
      // TODO: other codims

      virtual int subIndex0 (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().template subIndex<0>(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }

      virtual int subIndex1 (const typename GridImp::Traits::template Codim<1>::Entity& e, int i, int codim) const override
      {
        return impl().template subIndex<1>(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }

      virtual int subIndexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e, int i, int codim) const override
      {
        return impl().template subIndex<dim>(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }
      // TODO: other codims

      virtual std::size_t size (int codim) const override { return impl().size(codim); }
      virtual std::size_t size (GeometryType type) const override { return impl().size(type); }
      virtual Types types (int codim) const override { return impl().types(codim); }

      virtual bool contains0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().contains(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl()
        );
      }

      virtual bool contains1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const override
      {
        return impl().contains(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl()
        );
      }

      virtual bool containsDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const override
      {
        return impl().contains(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl()
        );
      }
      // TODO: other codims

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      const I& impl_;
    };
    // VIRTUALIZATION END


  public:
    template< class ImplLevelIndexSet >
    explicit VirtualizedGridLevelIndexSet(const ImplLevelIndexSet& implLevelIndexSet)
    : impl_( new Implementation<ImplLevelIndexSet>( implLevelIndexSet ) )
    {}

    VirtualizedGridLevelIndexSet(const VirtualizedGridLevelIndexSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLevelIndexSet ( VirtualizedGridLevelIndexSet && ) = default;

    VirtualizedGridLevelIndexSet& operator=(const VirtualizedGridLevelIndexSet& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      if constexpr (codim == 0)
        return impl_->index0(e);
      if constexpr (codim == 1)
        return impl_->index1(e);
      if constexpr (codim == dim)
        return impl_->indexDim(e);
    }


    //! get index of subEntity of a codim 0 entity
    template<int cc>
    int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      if constexpr (cc == 0)
        return impl_->subIndex0(e, i, codim);
      if constexpr (cc == 1)
        return impl_->subIndex1(e, i, codim);
      if constexpr (cc == dim)
        return impl_->subIndexDim(e, i, codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (int codim) const {
      return impl_->size(codim);
    }


    //! get number of entities of given codim, type and on this level
    std::size_t size (GeometryType type) const
    {
      return impl_->size(type);
    }


    /** \brief Deliver all geometry types used in this grid */
    const Types& geomTypes (int codim) const
    {
      return impl_->geomTypes(codim);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      return impl_->types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      static constexpr int cc = EntityType::codimension;
      if constexpr (cc == 0)
        return impl_->contains0(e);
      if constexpr (cc == 1)
        return impl_->contains1(e);
      if constexpr (cc == dim)
        return impl_->containsDim(e);
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
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual int index0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual int index1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const = 0;
      virtual int indexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const = 0;
      virtual int subIndex0 (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
      virtual int subIndex1 (const typename GridImp::Traits::template Codim<1>::Entity& e, int i, int codim) const = 0;
      virtual int subIndexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e, int i, int codim) const = 0;
      virtual std::size_t size (int codim) const = 0;
      virtual std::size_t size (GeometryType type) const = 0;
      virtual Types types (int codim) const = 0;
      virtual bool contains0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual bool contains1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const = 0;
      virtual bool containsDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename VirtualizedGridEntity<0, dim, GridImp>::template Implementation<typename I::template Codim<0>::Entity> ImplEntity0;
      typedef typename VirtualizedGridEntity<1, dim, GridImp>::template Implementation<typename I::template Codim<1>::Entity> ImplEntity1;
      typedef typename VirtualizedGridEntity<dim, dim, GridImp>::template Implementation<typename I::template Codim<dim>::Entity> ImplEntityDim;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual int index0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().index(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl()
        );
      }

      virtual int index1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const override
      {
        return impl().index(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl()
        );
      }

      virtual int indexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const override
      {
        return impl().template index(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl()
        );
      }
      // TODO: other codims

      virtual int subIndex0 (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().template subIndex<0>(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }

      virtual int subIndex1 (const typename GridImp::Traits::template Codim<1>::Entity& e, int i, int codim) const override
      {
        return impl().template subIndex<1>(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }

      virtual int subIndexDim (const typename GridImp::Traits::template Codim<dim>::Entity& e, int i, int codim) const override
      {
        return impl().template subIndex<dim>(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }
      // TODO: other codims

      virtual std::size_t size (int codim) const override { return impl().size(codim); }
      virtual std::size_t size (GeometryType type) const override { return impl().size(type); }
      virtual Types types (int codim) const override { return impl().types(codim); }

      virtual bool contains0 (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().contains(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl()
        );
      }

      virtual bool contains1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const override
      {
        return impl().contains(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl()
        );
      }

      virtual bool containsDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const override
      {
        return impl().contains(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl()
        );
      }
      // TODO: other codims

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      const I& impl_;
    };
    // VIRTUALIZATION END


  public:
    template< class ImplLeafIndexSet >
    explicit VirtualizedGridLeafIndexSet(const ImplLeafIndexSet& implLeafIndexSet)
    : impl_( new Implementation<ImplLeafIndexSet>( implLeafIndexSet ) )
    {}

    VirtualizedGridLeafIndexSet(const VirtualizedGridLeafIndexSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLeafIndexSet ( VirtualizedGridLeafIndexSet && ) = default;

    VirtualizedGridLeafIndexSet& operator=(const VirtualizedGridLeafIndexSet& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }


    //! get index of an entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int codim>
    int index (const typename GridImp::template Codim<codim>::Entity& e) const
    {
      if constexpr (codim == 0)
        return impl_->index0(e);
      if constexpr (codim == 1)
        return impl_->index1(e);
      if constexpr (codim == dim)
        return impl_->indexDim(e);
    }


    //! get index of subEntity of a codim 0 entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int cc>
    int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      if constexpr (cc == 0)
        return impl_->subIndex0(e, i, codim);
      if constexpr (cc == 1)
        return impl_->subIndex1(e, i, codim);
      if constexpr (cc == dim)
        return impl_->subIndexDim(e, i, codim);
    }


    //! get number of entities of given type
    std::size_t size (GeometryType type) const
    {
      return impl_->size(type);
    }


    //! get number of entities of given codim
    std::size_t size (int codim) const
    {
      return impl_->size(codim);
    }


    /** \brief Deliver all geometry types used in this grid */
    const Types& geomTypes (int codim) const
    {
      return impl_->geomTypes(codim);
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      return impl_->types(codim);
    }

    /** \brief Return true if the given entity is contained in the index set */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
      static constexpr int cc = EntityType::codimension;
      if constexpr (cc == 0)
        return impl_->contains0(e);
      if constexpr (cc == 1)
        return impl_->contains1(e);
      if constexpr (cc == dim)
        return impl_->containsDim(e);
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
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual IdType id1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const = 0;
      virtual IdType idDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const = 0;
      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename VirtualizedGridEntity<0, dim, GridImp>::template Implementation<typename I::template Codim<0>::Entity> ImplEntity0;
      typedef typename VirtualizedGridEntity<1, dim, GridImp>::template Implementation<typename I::template Codim<1>::Entity> ImplEntity1;
      typedef typename VirtualizedGridEntity<dim, dim, GridImp>::template Implementation<typename I::template Codim<dim>::Entity> ImplEntityDim;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().template id<0>(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl()
        );
      }

      virtual IdType id1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const override
      {
        return impl().template id<1>(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl()
        );
      }

      virtual IdType idDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const override
      {
        return impl().template id<dim>(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl()
        );
      }

      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().subId(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      const I& impl_;
    };
    // VIRTUALIZATION END


  public:
    template< class ImplGlobalIdSet >
    explicit VirtualizedGridGlobalIdSet(const ImplGlobalIdSet& implGlobalIdSet)
    : impl_( new Implementation<ImplGlobalIdSet>( implGlobalIdSet ) )
    {}

    VirtualizedGridGlobalIdSet(const VirtualizedGridGlobalIdSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridGlobalIdSet ( VirtualizedGridGlobalIdSet && ) = default;

    VirtualizedGridGlobalIdSet& operator=(const VirtualizedGridGlobalIdSet& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }


    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      if constexpr (cd == 0)
        return impl_->id(e);
      if constexpr (cd == 1)
        return impl_->id1(e);
      if constexpr (cd == dim)
        return impl_->idDim(e);
    }


    //! get id of subEntity
    IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      // Return sub id of the host entity
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
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual IdType id1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const = 0;
      virtual IdType idDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const = 0;
      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename VirtualizedGridEntity<0, dim, GridImp>::template Implementation<typename I::template Codim<0>::Entity> ImplEntity0;
      typedef typename VirtualizedGridEntity<1, dim, GridImp>::template Implementation<typename I::template Codim<1>::Entity> ImplEntity1;
      typedef typename VirtualizedGridEntity<dim, dim, GridImp>::template Implementation<typename I::template Codim<dim>::Entity> ImplEntityDim;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().template id<0>(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl()
        );
      }

      virtual IdType id1 (const typename GridImp::Traits::template Codim<1>::Entity& e) const override
      {
        return impl().template id<1>(
          dynamic_cast<const ImplEntity1*>(e.impl().impl_.get())->impl()
        );
      }

      virtual IdType idDim (const typename GridImp::Traits::template Codim<dim>::Entity& e) const override
      {
        return impl().template id<dim>(
          dynamic_cast<const ImplEntityDim*>(e.impl().impl_.get())->impl()
        );
      }

      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().subId(
          dynamic_cast<const ImplEntity0*>(e.impl().impl_.get())->impl(),
          i, codim
        );
      }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      const I& impl_;
    };
    // VIRTUALIZATION END


  public:
    template< class ImplLocalIdSet >
    explicit VirtualizedGridLocalIdSet(const ImplLocalIdSet& implLocalIdSet)
    : impl_( new Implementation<ImplLocalIdSet>( implLocalIdSet ) )
    {}

    VirtualizedGridLocalIdSet(const VirtualizedGridLocalIdSet& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridLocalIdSet ( VirtualizedGridLocalIdSet && ) = default;

    VirtualizedGridLocalIdSet& operator=(const VirtualizedGridLocalIdSet& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }


    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      if constexpr (cd == 0)
        return impl_->id(e);
      if constexpr (cd == 1)
        return impl_->id1(e);
      if constexpr (cd == dim)
        return impl_->idDim(e);
    }


    //! get id of subEntity
    IdType subId (const typename GridImp::template Codim<0>::Entity& e, int i, int codim) const
    {
      // Return sub id of the host entity
      return impl_->subId(e, i, codim);
    }


    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune


#endif
