// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRIDLEAFITERATOR_HH
#define DUNE_VIRTUALIZEDGRIDLEAFITERATOR_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/polymorphicgrid/common/iterator.hh>

/** \file
 * \brief The PolymorphicEntityIterator class
 */

namespace Dune {

  /** \brief Iterator over all entities of a given codimension and level of a grid.
   *  \ingroup Polymorphic
   */
  template<int codim, class GridImp>
  class PolymorphicEntityIterator
      : public Polymorphic::IteratorDefinition<
          PolymorphicEntityIterator<codim,GridImp>,
          PolymorphicEntity<codim,GridImp::dimension,GridImp>>::Base
  {
    using Self = PolymorphicEntityIterator;
    using EntityImp = PolymorphicEntity<codim,GridImp::dimension,GridImp>;
    using Definition = Polymorphic::IteratorDefinition<Self,EntityImp>;
    using Base = typename Definition::Base;

  public:
    enum {codimension = codim};

  public:
    using Entity = typename GridImp::template Codim<codim>::Entity;

    PolymorphicEntityIterator () = default;

    template <class Impl, disableCopyMove<PolymorphicEntityIterator,Impl> = 0>
    PolymorphicEntityIterator (Impl&& impl)
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
