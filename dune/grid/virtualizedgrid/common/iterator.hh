// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_POLYMORPHIC_COMMON_ITERATOR_HH
#define DUNE_POLYMORPHIC_COMMON_ITERATOR_HH

#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

/** \file
 * \brief The VirtualizedGridIntersectionIterator class
 */

namespace Dune::Polymorphic {

  template<class Iterator, class Reference>
  struct IteratorDefinition
  {
    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool equals (const Iterator&) const = 0;
      virtual void increment () = 0;
      virtual Reference dereference () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      using Wrapped = typename Wrapper::Wrapped;

      bool equals (const Iterator& other) const final
      {
        return this->get() == Polymorphic::asWrapped<Wrapped>(other);
      }

      void increment () final
      {
        ++(this->get());
      }

      Reference dereference () const final
      {
        return Reference{ *(this->get()) };
      }
    };

    using Base = TypeErasureBase<Interface, Implementation>;
  };

}  // namespace Dune::Polymorphic

#endif // DUNE_POLYMORPHIC_COMMON_ITERATOR_HH
