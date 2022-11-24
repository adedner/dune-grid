// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_ALBERTA_STRUCTUREDGRIDFACTORY_HH
#define DUNE_ALBERTA_STRUCTUREDGRIDFACTORY_HH

/** \file
 *  \author Simon Praetorius
 *  \brief  specialization of the generic StructuredGridFactory for AlbertaGrid
 */

#include <array>
#include <memory>
#include <vector>
#include <tuple>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/multiindex.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#if HAVE_ALBERTA

namespace Dune
{
  // Forward Declarations
  template <class Grid>
  class StructuredGridFactory;

  template <int dim, int dimworld>
  class AlbertaGrid;

  /** \brief specialization of the generic StructuredGridFactory for AlbertaGrid
   *
   * The simplex grid generated by the default StructuredGridFactory is not compatible
   * with the Alberta vertex numbering. This might lead to infinite recursions when
   * refining such grid. This specialization of the StructuredGridFactory solves this issue
   * by
   *
   * 1. generating a structured cube grid
   * 2. subdividing each cube into 2 (2d) or 6 (3d) simplices following a compatible numbering
   *    scheme for AlbertaGrid.
   *
   **/

  template <int dim, int dimworld>
  class StructuredGridFactory<AlbertaGrid<dim,dimworld>>
  {
    using GridType = AlbertaGrid<dim,dimworld>;

  protected:
    using ctype = typename GridType::ctype;

    // Insert new elements into the grid by splitting the given cube with corners `vertices` into
    // a corresponding number of simplices.
    static void insertElement (GridFactory<GridType>& factory,
      const GeometryType& type,
      const std::vector<unsigned int>& vertices)
    {
      // triangulation of reference cube
      static const auto reference_cubes = std::make_tuple(
        /*1d*/ std::array<std::array<int,2>,1>{
          std::array<int,2>{0,1}},
        /*2d*/ std::array<std::array<int,3>,2>{
          std::array<int,3>{3,0,1}, std::array<int,3>{0,3,2}},
        /*3d*/ std::array<std::array<int,4>,6>{
          std::array<int,4>{0,7,3,1}, std::array<int,4>{0,7,5,1},
          std::array<int,4>{0,7,5,4}, std::array<int,4>{0,7,3,2},
          std::array<int,4>{0,7,6,2}, std::array<int,4>{0,7,6,4}} );

      const auto& simplices = std::get<dim-1>(reference_cubes);
      std::vector<unsigned int> corners(dim+1);
      for (const auto& simplex : simplices) {
        for (std::size_t i = 0; i < simplex.size(); ++i)
          corners[i] = vertices[simplex[i]];

        factory.insertElement(type, corners);
      }
    }

    // Insert a structured set of vertices into the factory
    static void insertVertices (GridFactory<GridType>& factory,
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& vertices)
    {
      FactoryUtilities::MultiIndex<dim> index(vertices);

      // Compute the total number of vertices to be created
      int numVertices = index.cycle();

      // Create vertices
      for (int i=0; i<numVertices; i++, ++index) {

        // scale the multiindex to obtain a world position
        FieldVector<double,dimworld> pos(0);
        for (int j=0; j<dim; j++)
          pos[j] = lowerLeft[j] + index[j] * (upperRight[j]-lowerLeft[j])/(vertices[j]-1);
        for (int j=dim; j<dimworld; j++)
          pos[j] = lowerLeft[j];

        factory.insertVertex(pos);

      }

    }

    // Compute the index offsets needed to move to the adjacent vertices
    // in the different coordinate directions
    static std::array<unsigned int, dim> computeUnitOffsets (
      const std::array<unsigned int,dim>& vertices)
    {
      std::array<unsigned int, dim> unitOffsets;
      if (dim>0)        // paranoia
        unitOffsets[0] = 1;

      for (int i=1; i<dim; i++)
        unitOffsets[i] = unitOffsets[i-1] * vertices[i-1];

      return unitOffsets;
    }

  public:
    /** \brief Create a structured simplex grid for AlbertaGrid
     *
     * This works in dimensions 1-3. Each cube is split into into dim! simplices by constructing
     * macro triangulation, see \ref insertElement for the local simplex numbering.
     *
     * \param factory     Grid factory used for creating the grid
     * \param lowerLeft   Lower left corner of the grid
     * \param upperRight  Upper right corner of the grid
     * \param elements    Number of elements in each coordinate direction
     **/
    static void createSimplexGrid (GridFactory<GridType>& factory,
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      if (factory.comm().rank() == 0)
      {
        // Insert uniformly spaced vertices
        std::array<unsigned int,dim> vertices = elements;
        for (std::size_t i=0; i<vertices.size(); ++i)
          vertices[i]++;

        insertVertices(factory, lowerLeft, upperRight, vertices);

        // Compute the index offsets needed to move to the adjacent
        // vertices in the different coordinate directions
        std::array<unsigned int, dim> unitOffsets =
          computeUnitOffsets(vertices);

        // Compute an element template (the cube at (0,...,0).  All
        // other cubes are constructed by moving this template around
        unsigned int nCorners = 1<<dim;

        std::vector<unsigned int> cornersTemplate(nCorners,0);
        for (std::size_t i=0; i<nCorners; ++i)
          for (int j=0; j<dim; ++j)
            if ( i & (1<<j) )
              cornersTemplate[i] += unitOffsets[j];

        // Insert elements
        FactoryUtilities::MultiIndex<dim> index(elements);

        // Compute the total number of elements to be created
        int numElements = index.cycle();

        for (int i=0; i<numElements; ++i, ++index) {

          // 'base' is the index of the lower left element corner
          unsigned int base = 0;
          for (int j=0; j<dim; j++)
            base += index[j] * unitOffsets[j];

          // insert new element
          std::vector<unsigned int> corners = cornersTemplate;
          for (std::size_t j=0; j<corners.size(); ++j)
            corners[j] += base;

          insertElement(factory, GeometryTypes::simplex(dim), corners);
        }
      }
    }

    /** \brief Create a structured simplex grid for AlbertaGrid
     *
     * This works in dimensions 1-3. Each cube is split into into dim! simplices by constructing
     * macro triangulation, see \ref insertElement for the local simplex numbering.
     *
     * \param lowerLeft   Lower left corner of the grid
     * \param upperRight  Upper right corner of the grid
     * \param elements    Number of elements in each coordinate direction
     **/
    static std::unique_ptr<GridType> createSimplexGrid (
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      GridFactory<GridType> factory;
      createSimplexGrid(factory, lowerLeft, upperRight, elements);
      return std::unique_ptr<GridType>(factory.createGrid());
    }


    static void createCubeGrid (GridFactory<GridType>& factory,
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      DUNE_THROW(Dune::NotImplemented,
        "Cube grids are not supported by AlbertaGrid. Use createSimplexGrid instead.");
    }

    static std::unique_ptr<GridType> createCubeGrid (
      const FieldVector<ctype,dimworld>& lowerLeft,
      const FieldVector<ctype,dimworld>& upperRight,
      const std::array<unsigned int,dim>& elements)
    {
      DUNE_THROW(Dune::NotImplemented,
        "Cube grids are not supported by AlbertaGrid. Use createSimplexGrid instead.");
      return nullptr;
    }
  };

} // end namespace Dune

#endif // HAVE_ALBERTA

#endif // DUNE_ALBERTA_STRUCTUREDGRIDFACTORY_HH
