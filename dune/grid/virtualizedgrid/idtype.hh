// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZED_GRID_IDTYPE_HH
#define DUNE_VIRTUALIZED_GRID_IDTYPE_HH

/**
 * \file
 * \brief The VirtualizedGridIdType class
 */


namespace Dune {


  /**
   * \brief The IdType class provides a virtualized id type.
   * \ingroup VirtualizedGrid
   *
   */
  class VirtualizedGridIdType
  {

  public:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual bool operator== (const VirtualizedGridIdType& other) const = 0;
      virtual bool operator!= (const VirtualizedGridIdType& other) const = 0;
      virtual bool operator< (const VirtualizedGridIdType& other) const = 0;
      virtual std::string str() const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::move(i) ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual bool operator== (const VirtualizedGridIdType& other) const override
      {
        return impl() == dynamic_cast<const Implementation<I>&>(*other.impl_).impl();
      }

      virtual bool operator!= (const VirtualizedGridIdType& other) const override
      {
        return !operator==(other);
      }

      virtual bool operator< (const VirtualizedGridIdType& other) const override
      {
        return impl() < dynamic_cast<const Implementation<I>&>(*other.impl_).impl();
      }

      virtual std::string str() const override
      {
        std::stringstream ss;
        ss << impl() << std::endl;
        return ss.str();
      }

      const auto &impl () const { return impl_; }
    private:
      auto &impl () { return impl_; }

      const I impl_;
    };
    // VIRTUALIZATION END

    VirtualizedGridIdType()
    {}

    template< class ImplIdType >
    VirtualizedGridIdType(ImplIdType&& implIdType)
    : impl_( new Implementation<ImplIdType>( std::move(implIdType) ) )
    {}

    VirtualizedGridIdType(const VirtualizedGridIdType& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridIdType ( VirtualizedGridIdType && ) = default;

    VirtualizedGridIdType& operator=(const VirtualizedGridIdType& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      return *this;
    }

    VirtualizedGridIdType& operator=( VirtualizedGridIdType&& ) = default;

    bool operator==(const VirtualizedGridIdType& other) const
    {
      return impl_->operator==(other);
    }

    bool operator!=(const VirtualizedGridIdType& other) const
    {
      return impl_->operator!=(other);
    }

    bool operator<(const VirtualizedGridIdType& other) const
    {
      return impl_->operator<(other);
    }

    std::string str() const
    {
      return impl_->str();
    }

    std::unique_ptr<Interface> impl_;
  };

  inline std::ostream &operator<< ( std::ostream &out, const VirtualizedGridIdType &idtype )
  {
    return out << idtype.str();
  }

} // namespace Dune

#endif  // #define DUNE_VIRTUALIZED_GRID_IDTYPE_HH
