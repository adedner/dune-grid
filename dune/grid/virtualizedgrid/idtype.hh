// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZED_GRID_IDTYPE_HH
#define DUNE_VIRTUALIZED_GRID_IDTYPE_HH

/**
 * \file
 * \brief The VirtualizedGridIdType class
 */

#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

namespace Dune {

  // forward declaration
  class VirtualizedGridIdType;

  struct VirtualizedGridIdTypeDefinition
  {
    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool operator== (const VirtualizedGridIdType&) const = 0;
      virtual bool operator!= (const VirtualizedGridIdType&) const = 0;
      virtual bool operator< (const VirtualizedGridIdType&) const = 0;
      virtual bool operator<= (const VirtualizedGridIdType&) const = 0;
      virtual std::size_t hash () const = 0;
      virtual std::string str () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      using Wrapped = typename Wrapper::Wrapped;

      bool operator== (const VirtualizedGridIdType& other) const final
      {
        return this->get() == Polymorphic::asWrapped<Wrapped>(other);
      }

      bool operator!= (const VirtualizedGridIdType& other) const final
      {
        return this->get() != Polymorphic::asWrapped<Wrapped>(other);
      }

      bool operator< (const VirtualizedGridIdType& other) const final
      {
        return this->get() < Polymorphic::asWrapped<Wrapped>(other);
      }

      bool operator<= (const VirtualizedGridIdType& other) const final
      {
        return this->get() <= Polymorphic::asWrapped<Wrapped>(other);
      }

      std::string str () const final
      {
        std::stringstream ss;
        ss << this->get() << std::endl;
        return ss.str();
      }

      std::size_t hash () const final
      {
        return std::hash<Wrapped>()(this->get());
      }
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };

  /**
   * \brief The IdType class provides a virtualized id type.
   * \ingroup VirtualizedGrid
   *
   */
  class VirtualizedGridIdType :
      public VirtualizedGridIdTypeDefinition::Base
  {
    using Definition = VirtualizedGridIdTypeDefinition;
    using Base = typename Definition::Base;

  public:
    VirtualizedGridIdType () = default;

    /**
     * \brief Create IdType from implementation
     */
    template <class Impl, disableCopyMove<VirtualizedGridIdType,Impl> = 0>
    VirtualizedGridIdType (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    bool operator== (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator==(other);
    }

    bool operator!= (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator!=(other);
    }

    bool operator< (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator<(other);
    }

    bool operator<= (const VirtualizedGridIdType& other) const
    {
      return this->asInterface().operator<=(other);
    }

    std::string str () const
    {
      return this->asInterface().str();
    }

    std::size_t hash () const
    {
      return this->asInterface().hash();
    }
  };

  inline std::ostream& operator<< (std::ostream &out, const VirtualizedGridIdType& idtype)
  {
    return out << idtype.str();
  }

} // namespace Dune


namespace std
{
  template <> struct hash<Dune::VirtualizedGridIdType>
  {
    size_t operator() (const Dune::VirtualizedGridIdType& x) const
    {
      return x.hash();
    }
  };
}

#endif  // #define DUNE_VIRTUALIZED_GRID_IDTYPE_HH
