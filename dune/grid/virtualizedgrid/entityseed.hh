// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH
#define DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH

/**
 * \file
 * \brief The VirtualizedGridEntitySeed class
 */

#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

namespace Dune {

  struct VirtualizedGridEntitySeedDefinition
  {
    struct Interface
    {
      virtual ~Interface () = default;
      virtual bool isValid () const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;
      bool isValid () const final { return this->get().isValid(); }
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };


  /**
   * \brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup VirtualizedGrid
   *
   */
  template<int codim, class GridImp>
  class VirtualizedGridEntitySeed
      : public VirtualizedGridEntitySeedDefinition::Base
  {
    using Definition = VirtualizedGridEntitySeedDefinition;
    using Base = typename Definition::Base;

  protected:

    // Entity type of the grid
    using Entity = typename GridImp::Traits::template Codim<codim>::Entity;

  public:

    enum {codimension = codim};

    /**
     * \brief Construct an empty (i.e. isValid() == false) seed.
     */
    VirtualizedGridEntitySeed () = default;

    /**
     * \brief Create EntitySeed from implementation entity
     */
    template <class Impl, disableCopyMove<VirtualizedGridEntitySeed,Impl> = 0>
    VirtualizedGridEntitySeed (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    /**
     * \brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid () const
    {
      return this->wrapped_ && this->asInterface().isValid();
    }
  };

} // namespace Dune

#endif  // #define DUNE_VIRTUALIZED_GRID_ENTITY_SEED_HH
