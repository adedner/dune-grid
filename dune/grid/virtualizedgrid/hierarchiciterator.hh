// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRIDHIERITERATOR_HH
#define DUNE_VIRTUALIZEDGRIDHIERITERATOR_HH

/** \file
 * \brief The VirtualizedGridHierarchicIterator class
 */

#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

namespace Dune {

  template<class GridImp>
  struct VirtualizedGridHierarchicIteratorDefinition
  {
    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool equals (const VirtualizedGridHierarchicIterator<GridImp>&) const = 0;
      virtual void increment () = 0;
      virtual Entity dereference () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      bool equals (const VirtualizedGridHierarchicIterator<GridImp>& other) const final
      {
        return this->get() == Polymorphic::asWrapped<Wrapped>(other);
      }
      void increment () final { this->get().increment(); };
      Entity dereference () const final { return this->get().dereference(); };
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };

  //**********************************************************************
  //
  /** \brief Iterator over the descendants of an entity.
   * \ingroup VirtualizedGrid
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HierarchicIterator,
     starting from a given entity.
   */
  template<class GridImp>
  class VirtualizedGridHierarchicIterator :
      public VirtualizedGridHierarchicIteratorDefinition<GridImp>::Base
  {
    using Definition = VirtualizedGridHierarchicIteratorDefinition<GridImp>;
    using Base = typename Definition::Base;

  public:
    enum { codimension = 0 };

    using Entity = typename GridImp::template Codim<0>::Entity;

  public:
    /**
     * \brief Create HierarchicIterator from implementation
     */
    template <class Impl, disableCopyMove<VirtualizedGridHierarchicIterator,Impl> = 0>
    VirtualizedGridHierarchicIterator (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    void increment ()
    {
      this->asInterface().increment();
    }

    Entity dereference () const
    {
      return this->asInterface().dereference();
    }

    bool equals (const VirtualizedGridHierarchicIterator& other) const
    {
      return this->asInterface().equals(other);
    }
  };

}  // end namespace Dune

#endif
