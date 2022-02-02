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
    public IndexSet<GridImp, VirtualizedGridLevelIndexSet<GridImp>, int, std::vector< GeometryType >>
  {
  public:

    typedef std::vector< GeometryType > Types;

    enum {dim = GridImp::dimension};

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual int index (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
      virtual std::size_t size (int codim) const = 0;
      virtual std::size_t size (GeometryType type) const = 0;
      virtual Types types (int codim) const = 0;
      virtual bool contains (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename I::template Codim<0>::Entity ImplEntity;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual int index (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().index(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_))
        );
      } // TODO: other codims

      virtual int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().subIndex(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_)), i, codim
        );
      }  // TODO: other codims

      virtual std::size_t size (int codim) const override { return impl().size(codim); }
      virtual std::size_t size (GeometryType type) const override { return impl().size(type); }
      virtual Types types (int codim) const override { return impl().types(codim); }

      virtual bool contains (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().contains(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_))
        );
      } // TODO: other codims

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
    }

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      return impl_->template index<codim>(e);
    }


    //! get index of subEntity of a codim 0 entity
    template<int cc>
    int subIndex (const typename GridImp::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      return impl_->subIndex(e, i, codim);
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
      return impl_->contains(e);
    }

    std::unique_ptr<Interface> impl_;
  };


  template<class GridImp>
  class VirtualizedGridLeafIndexSet :
    public IndexSet<GridImp, VirtualizedGridLeafIndexSet<GridImp>, int, std::vector< GeometryType >>
  {

  public:
    typedef std::vector< GeometryType > Types;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    enum {dim = std::remove_const<GridImp>::type::dimension};

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual int index (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
      virtual int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
      virtual std::size_t size (int codim) const = 0;
      virtual std::size_t size (GeometryType type) const = 0;
      virtual Types types (int codim) const = 0;
      virtual bool contains (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename I::template Codim<0>::Entity ImplEntity;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual int index (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().index(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_))
        );
      } // TODO: other codims

      virtual int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().subIndex(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_)), i, codim
        );
      }  // TODO: other codims

      virtual std::size_t size (int codim) const override { return impl().size(codim); }
      virtual std::size_t size (GeometryType type) const override { return impl().size(type); }
      virtual Types types (int codim) const override { return impl().types(codim); }
      virtual bool contains (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().contains(*dynamic_cast<ImplEntity*>(&(*e.impl().impl_)));
      } // TODO: other codims

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
    }


    //! get index of an entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int codim>
    int index (const typename std::remove_const<GridImp>::type::template Codim<codim>::Entity& e) const
    {
      return impl_->template index<codim>(e);
    }


    //! get index of subEntity of a codim 0 entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    template<int cc>
    int subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      return impl_->subIndex(e,i, codim);
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
      return impl_->contains(e);
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

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0; // TODO: other codims
      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename I::template Codim<0>::Entity ImplEntity;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().template id<0>(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_))
        );
      } // TODO: other codims
      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().subId(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_)), i, codim
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
    }


    //! get id of an entity
    template<int cd>
    IdType id (const typename GridImp::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return impl_->id(e);
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

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const = 0; // TODO: other codims
      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      typedef typename I::template Codim<0>::Entity ImplEntity;

      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }
      virtual IdType id (const typename GridImp::Traits::template Codim<0>::Entity& e) const override
      {
        return impl().template id<0>(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_))
        );
      } // TODO: other codims
      virtual IdType subId (const typename GridImp::Traits::template Codim<0>::Entity& e, int i, int codim) const override
      {
        return impl().subId(
          *dynamic_cast<ImplEntity*>(&(*e.impl().impl_)), i, codim
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
    }


    //! get id of an entity
    template<int cd>
    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return impl_->id(e.impl().impl_);
    }


    //! get id of subEntity
    IdType subId (const typename std::remove_const<GridImp>::type::template Codim<0>::Entity& e, int i, int codim) const
    {
      // Return sub id of the host entity
      return impl_->subId(e.impl().impl_, i, codim);
    }


    std::unique_ptr<Interface> impl_;
  };


}  // namespace Dune


#endif
