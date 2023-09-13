// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZED_GRID_IDTYPE_HH
#define DUNE_VIRTUALIZED_GRID_IDTYPE_HH

/**
 * \file
 * \brief The PolymorphicIdType class
 */

#include <dune/grid/polymorphicgrid/common/typeerasure.hh>

namespace Dune {

  // forward declaration
  class PolymorphicIdType;

  struct PolymorphicIdTypeDefinition
  {
    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool operator== (const PolymorphicIdType&) const = 0;
      virtual bool operator!= (const PolymorphicIdType&) const = 0;
      virtual bool operator< (const PolymorphicIdType&) const = 0;
      virtual bool operator<= (const PolymorphicIdType&) const = 0;
      virtual std::size_t hash () const = 0;
      virtual std::string str () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      using Wrapped = typename Wrapper::Wrapped;

      bool operator== (const PolymorphicIdType& other) const final
      {
        return this->get() == Polymorphic::asWrapped<Wrapped>(other);
      }

      bool operator!= (const PolymorphicIdType& other) const final
      {
        return this->get() != Polymorphic::asWrapped<Wrapped>(other);
      }

      bool operator< (const PolymorphicIdType& other) const final
      {
        return this->get() < Polymorphic::asWrapped<Wrapped>(other);
      }

      bool operator<= (const PolymorphicIdType& other) const final
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
   * \brief The IdType class provides a polymorphic id type.
   * \ingroup Polymorphic
   *
   */
  class PolymorphicIdType :
      public PolymorphicIdTypeDefinition::Base
  {
    using Definition = PolymorphicIdTypeDefinition;
    using Base = typename Definition::Base;

  public:
    PolymorphicIdType () = default;

    /**
     * \brief Create IdType from implementation
     */
    template <class Impl, disableCopyMove<PolymorphicIdType,Impl> = 0>
    PolymorphicIdType (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    bool operator== (const PolymorphicIdType& other) const
    {
      return this->asInterface().operator==(other);
    }

    bool operator!= (const PolymorphicIdType& other) const
    {
      return this->asInterface().operator!=(other);
    }

    bool operator< (const PolymorphicIdType& other) const
    {
      return this->asInterface().operator<(other);
    }

    bool operator<= (const PolymorphicIdType& other) const
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

  inline std::ostream& operator<< (std::ostream &out, const PolymorphicIdType& idtype)
  {
    return out << idtype.str();
  }

} // namespace Dune


namespace std
{
  template <> struct hash<Dune::PolymorphicIdType>
  {
    size_t operator() (const Dune::PolymorphicIdType& x) const
    {
      return x.hash();
    }
  };
}

#endif  // #define DUNE_VIRTUALIZED_GRID_IDTYPE_HH
