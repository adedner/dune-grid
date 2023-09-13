// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DGFPARSERVIRTUALIZED_HH
#define DUNE_DGFPARSERVIRTUALIZED_HH

#include <dune/grid/polymorphicgrid.hh>
#include "dgfparser.hh"

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class GridImp, class IntersectionImp >
  class Intersection;


  /*!
   * \brief Dummy grid factory for Polymorphic
   */
  template <int dim, int dimw>
  struct DGFGridFactory< PolymorphicGrid<dim, dimw> >
  {
    typedef PolymorphicGrid<dim, dimw> Grid;
    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;

  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for Polymorphic is not available!");
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for Polymorphic is not available!");
    }

    Grid *grid() const
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for Polymorphic is not available!");
      return nullptr;
    }

  };

  template <int dim, int dimw>
  struct DGFGridInfo< PolymorphicGrid<dim, dimw> >
  {
    static int refineStepsForHalf()
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for Polymorphic is not available!");
      return 1;
    }

    static double refineWeight()
    {
      DUNE_THROW(NotImplemented, "DGFGridFactory for Polymorphic is not available!");
      return 1;
    }
  };

}
#endif // #ifndef DUNE_DGFPARSERVIRTUALIZED_HH
