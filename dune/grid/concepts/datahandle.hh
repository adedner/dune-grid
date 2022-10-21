// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_DATAHANDLE_HH
#define DUNE_GRID_CONCEPTS_DATAHANDLE_HH

#include <dune/grid/concepts/archetypes/datahandle.hh>
#include <dune/grid/concepts/archetypes/entity.hh>

namespace Dune::Concept {

template <class MB, class DataType>
concept MessageBuffer = requires(MB buffer, DataType data)
{
  write(data);
  read(data);
};

static_assert(MessageBuffer< Archetypes::MessageBuffer<double>, double >);


template <class DH>
concept CommDataHandle = requires(DH handle, const Archetypes::Entity<2,0>& entity)
{
  typename DH::DataType;

  { handle.contains(/*dim*/ 0, /*codim*/ 0) } -> std::convertible_to<bool>;
  { handle.fixedSize(/*dim*/ 0, /*codim*/ 0) } -> std::convertible_to<bool>;
  { handle.size(entity) } -> std::integral;

  requires requires(Archetypes::MessageBuffer<typename DH::DataType> buffer)
  {
    handle.gather(buffer, entity);
    handle.scatter(buffer, entity, /*size*/ 0u);
  }
};

static_assert(CommDataHandle< Archetypes::CommDataHandle<double> >);

} // end namespace Dune::Concept


#endif // DUNE_GRID_CONCEPTS_DATAHANDLE_HH
