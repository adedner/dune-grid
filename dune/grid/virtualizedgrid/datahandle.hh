// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_VERTUALIZEDGRID_DATAHANDLE_HH
#define DUNE_GRID_VERTUALIZEDGRID_DATAHANDLE_HH

#include <utility>

#include <dune/common/visibility.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/virtualizedgrid/entity.hh>

namespace Dune
{
  struct VirtualizedMessageBuffer
  {
    std::function<void(const std::byte*)> writeImpl;
    std::function<void(std::byte*)> readImpl;

    template<class T>
    void write (const T& data)
    {
      const std::byte* byte_data = reinterpret_cast<const std::byte*>(&data);
      writeImpl(byte_data);
    }

    template<class T>
    void read (T& data)
    {
      std::byte* byte_data = reinterpret_cast<std::byte*>(&data);
      readImpl(byte_data);
    }
  };

  template <class Data, class GridImp>
  class VirtualizedCommDataHandle
  {
  public:
    using DataType = Data;
    using MessageBuffer = VirtualizedMessageBuffer;

  private:
    // VIRTUALIZATION BEGIN
    template<int codim>
    struct InterfaceCodim
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;

      virtual ~InterfaceCodim () = default;
      virtual std::size_t size (Codim<codim>, const Entity& e) const = 0;
      virtual void gather (Codim<codim>, MessageBuffer& buff, const Entity& e) const = 0;
      virtual void scatter (Codim<codim>, MessageBuffer& buff, const Entity& e, std::size_t n) = 0;
    };

    template<class Derived, class I, int codim>
    struct DUNE_PRIVATE ImplementationCodim
        : virtual InterfaceCodim<codim>
    {
      using Entity = typename GridImp::Traits::template Codim<codim>::Entity;
      using EntityImpl = typename VirtualizedGridEntity<codim,GridImp::dimension,GridImp>::template Implementation<const typename std::decay_t<I>::template Codim<codim>::Entity>;

      std::size_t size (Codim<codim>, const Entity& e) const final {
        return derived().impl().size(upcast<EntityImpl>(e));
      }

      void gather (Codim<codim>, MessageBuffer& buff, const Entity& e) const final {
        derived().impl().gather(buff, upcast<EntityImpl>(e));
      }

      void scatter (Codim<codim>, MessageBuffer& buff, const Entity& e, std::size_t n) final {
        derived().impl().scatter(buff, upcast<EntityImpl>(e), n);
      }

    private:
      const Derived& derived () const { return static_cast<const Derived&>(*this); }
    };


    template<int... codims>
    struct InterfaceImpl
        : virtual InterfaceCodim<codims>...
    {
      virtual ~InterfaceImpl () = default;
      virtual InterfaceImpl *clone () const = 0;

      virtual bool contains (int dim, int codim) const = 0;
      virtual bool fixedSize (int dim, int codim) const = 0;

      using InterfaceCodim<codims>::size...;
      using InterfaceCodim<codims>::gather...;
      using InterfaceCodim<codims>::scatter...;
    };

    template<class I, int... codims>
    struct DUNE_PRIVATE ImplementationImpl final
        : virtual InterfaceImpl<codims...>
        , public ImplementationCodim<ImplementationImpl<I,codims...>, I, codims>...
    {
      ImplementationImpl (I&& i) : impl_(std::forward<I>(i)) {}
      ImplementationImpl* clone() const override { return new ImplementationImpl(*this); }

      bool contains (int dim, int codim) const override {
        return impl().contains(dim, codim);
      }

      bool fixedSize (int dim, int codim) const override {
        return impl().fixedSize(dim, codim);
      }

      const auto& impl () const { return impl_; }
      auto& impl () { return impl_; }

    private:
      I impl_;
    };
    // VIRTUALIZATION END

    template<class Seq> struct Interface_t;
    template<int... codims>
    struct Interface_t<std::integer_sequence<int,codims...>> {
      using type = InterfaceImpl<codims...>;
    };
    using Interface = typename Interface_t<std::make_integer_sequence<int,GridImp::dimension+1>>::type;

    template<class I, class Seq> struct Implementation_t;
    template<class I, int... codims>
    struct Implementation_t<I,std::integer_sequence<int,codims...>> {
      using type = ImplementationImpl<I,codims...>;
    };
    template<class I>
    using Implementation = typename Implementation_t<I,std::make_integer_sequence<int,GridImp::dimension+1>>::type;

  public:
    template<class Impl, class MessageBufferCreator>
    explicit VirtualizedCommDataHandle(Impl&& impl, const MessageBufferCreator& creator)
    : impl_(new Implementation<Impl>(std::forward<Impl>(impl)))
    , messageBufferCreator_(creator)
    {}

    VirtualizedCommDataHandle (const VirtualizedCommDataHandle& other)
    : impl_(other.impl_ ? other.impl_->clone() : nullptr)
    , messageBufferCreator_(other.messageBufferCreator_)
    {}

    VirtualizedCommDataHandle& operator= (const VirtualizedCommDataHandle& other) {
      impl_.reset(other.impl_ ? other.impl_->clone() : nullptr);
      messageBufferCreator_ = other.messageBufferCreator_;
      return *this;
    }

    VirtualizedCommDataHandle (VirtualizedCommDataHandle&&) = default;
    VirtualizedCommDataHandle& operator= (VirtualizedCommDataHandle&&) = default;


    bool contains (int dim, int codim) const {
      return impl_.contains(dim, codim);
    }

    bool fixedSize (int dim, int codim) const {
      return impl_.fixedSize(dim, codim);
    }

    template<class Entity>
    std::size_t size (const Entity& e) const {
      return impl_.size(Codim<Entity::codimension>{}, e);
    }

    template<class MessageBuffer, class Entity>
    void gather (MessageBuffer& buff, const Entity& e) const {
      VirtualizedMessageBuffer mb{
        [&buff](const std::byte* data) { buff.write(*reinterpret_cast<const Data*>(data)); },
        [&buff](std::byte* data) { buff.read(*reinterpret_cast<Data*>(data)); }
      };
      impl_.gather(Codim<Entity::codimension>{}, mb, e);
    }

    template<class MessageBuffer, class Entity>
    void scatter (MessageBuffer& buff, const Entity& e, std::size_t n) {
      VirtualizedMessageBuffer mb{
        [&buff](const std::byte* data) { buff.write(*reinterpret_cast<const Data*>(data)); },
        [&buff](std::byte* data) { buff.read(*reinterpret_cast<Data*>(data)); }
      };
      impl_.scatter(Codim<Entity::codimension>{}, mb, e, n);
    }

  private:
    std::unique_ptr<Interface> impl_;
  };

} // end namespace Dune

#endif // DUNE_GRID_VERTUALIZEDGRID_DATAHANDLE_HH
