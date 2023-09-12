// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRIDLEAFITERATOR_HH
#define DUNE_VIRTUALIZEDGRIDLEAFITERATOR_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/virtualizedgrid/common/iterator.hh>

/** \file
 * \brief The VirtualizedGridEntityIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup VirtualizedGrid
   */
  template<int codim, class GridImp>
  class VirtualizedGridEntityIterator
      : public Polymorphic::IteratorDefinition<
          VirtualizedGridEntityIterator<codim,GridImp>,
          VirtualizedGridEntity<codim,GridImp::dimension,GridImp>>::Base
  {
    using Self = VirtualizedGridEntityIterator;
    using EntityImp = VirtualizedGridEntity<codim,GridImp::dimension,GridImp>;
    using Definition = Polymorphic::IteratorDefinition<Self,EntityImp>;
    using Base = typename Definition::Base;

  public:
    enum {codimension = codim};

  public:
    using Entity = typename GridImp::template Codim<codim>::Entity;

    VirtualizedGridEntityIterator () = default;

    template <class Impl, disableCopyMove<VirtualizedGridEntityIterator,Impl> = 0>
    VirtualizedGridEntityIterator (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    //! prefix increment
    void increment ()
    {
      this->asInterface().increment();
    }

    //! dereferencing
    Entity dereference () const
    {
      return Entity{ this->asInterface().dereference() };
    }

    //! equality
    bool equals (const Self& other) const
    {
      return this->asInterface().equals(other);
    }
  };


}  // namespace Dune

#endif
