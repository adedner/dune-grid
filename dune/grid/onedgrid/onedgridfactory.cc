// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <dune/grid/onedgrid/onedgridfactory.hh>
#include <dune/grid/onedgrid/onedgridindexsets.hh>

namespace Dune {

template <int dimw, class ct>
GridFactory<OneDEmbeddedGrid<dimw,ct>>::
GridFactory() :
  factoryOwnsGrid_(true),
  vertexIndex_(0)
{
  grid_ = new OneDEmbeddedGrid<dimw,ct>;

  createBegin();
}

template <int dimw, class ct>
GridFactory<OneDEmbeddedGrid<dimw,ct>>::
GridFactory(Grid* grid) :
  factoryOwnsGrid_(false),
  vertexIndex_(0)
{
  grid_ = grid;

  createBegin();
}

template <int dimw, class ct>
GridFactory<OneDEmbeddedGrid<dimw,ct>>::
~GridFactory()
{
  if (grid_ && factoryOwnsGrid_)
    delete grid_;
}

template <int dimw, class ct>
void GridFactory<OneDEmbeddedGrid<dimw,ct>>::
insertVertex(const Dune::FieldVector<ct,dimw>& pos)
{
  vertexPositions_.insert(std::make_pair(pos, vertexIndex_++));
}

template <int dimw, class ct>
void GridFactory<OneDEmbeddedGrid<dimw,ct>>::
insertElement(const GeometryType& type,
              const std::vector<unsigned int>& vertices)
{
  if (type.dim() != 1)
    DUNE_THROW(GridError, "You cannot insert a " << type << " into a OneDGrid!");

  if (vertices.size() != 2)
    DUNE_THROW(GridError, "You cannot insert an element with " << vertices.size() << " vertices into a OneDGrid!");

  elements_.push_back(std::array<unsigned int,2>());
  elements_.back()[0] = vertices[0];
  elements_.back()[1] = vertices[1];
}

template <int dimw, class ct>
void GridFactory<OneDEmbeddedGrid<dimw,ct>>::
insertBoundarySegment(const std::vector<unsigned int>& vertices)
{
  if (vertices.size() != 1)
    DUNE_THROW(GridError, "OneDGrid BoundarySegments must have exactly one vertex.");

  boundarySegments_.push_back(vertices[0]);
}

template <int dimw, class ct>
void GridFactory<OneDEmbeddedGrid<dimw,ct>>::
insertBoundarySegment(const std::vector<unsigned int>& vertices,
                      const std::shared_ptr<BoundarySegment<1> > & /* boundarySegment */)
{
  insertBoundarySegment(vertices);
}

template <int dimw, class ct>
bool GridFactory<OneDEmbeddedGrid<dimw,ct>>::
wasInserted(const typename Grid::LeafIntersection& intersection) const
{
  bool inserted(false);
  const auto vtx(intersection.geometry().center());
  for(const auto& idx : boundarySegments_)
    if((vertexPositionsByIndex_[idx]-vtx).two_norm()<1.e-12)
    {
      inserted=true;
      break;
    }
  return inserted;
}

template <int dimw, class ct>
unsigned int GridFactory<OneDEmbeddedGrid<dimw,ct>>::
insertionIndex(const typename Grid::LeafIntersection& intersection) const
{
  unsigned int insertionIdx(0);
  const auto vtx(intersection.geometry().center());
  for(const auto& idx : boundarySegments_)
    if((vertexPositionsByIndex_[idx]-vtx).two_norm()<1.e-12)
      break;
    else
      ++insertionIdx;
  return insertionIdx;
}

template <int dimw, class ct>
std::unique_ptr<OneDEmbeddedGrid<dimw,ct>> GridFactory<OneDEmbeddedGrid<dimw,ct>>::
createGrid()
{
  // Prevent a crash when this method is called twice in a row
  // You never know who may do this...
  if (grid_==nullptr)
    return nullptr;

  // Assert that vertices are given
  assert (vertexPositions_.size() > 0);

  // Insert the vertices into the grid
  grid_->entityImps_.resize(1);
  for (const auto& vtx : vertexPositions_)
  {
    OneDEntityImp<0,dimw,ct> newVertex(0, vtx.first, grid_->getNextFreeId());

    newVertex.leafIndex_ = vtx.second;
    newVertex.levelIndex_ = vtx.second;

    grid_->vertices(0).push_back(newVertex);
  }

  // Fill the vector with the vertex positions accessible by index
  vertexPositionsByIndex_.resize(vertexPositions_.size());
  for (const auto& vtx : vertexPositions_)
    vertexPositionsByIndex_[vtx.second] = vtx.first;

  // Set the numbering of the boundary segments
  if (boundarySegments_.size() > 2)
    DUNE_THROW(GridError, "You cannot provide more than two boundary segments to a OneDGrid (it must be connected).");

  Impl::CompareFieldVector cmp;
  if (boundarySegments_.size() > 1
      && cmp(vertexPositions_.begin()->first, vertexPositionsByIndex_[boundarySegments_[0]]))
    grid_->reversedBoundarySegmentNumbering_ = true;

  // ///////////////////////////////////////////////////////////////////
  //   Insert the elements into the grid
  //
  // This is a 1d grid and it has to be connected. Hence we actually
  // know where the elements are, even without being told explicitly.
  // The only thing of interest are the indices.
  // ///////////////////////////////////////////////////////////////////

  // First sort elements by increasing position. That is how they are expected in the grid data structure
  std::map<Dune::FieldVector<ct,dimw>, std::pair<std::array<unsigned int, 2>, unsigned int>, Impl::CompareFieldVector > elementsByPosition;
  for (std::size_t i=0; i<elements_.size(); i++)
    elementsByPosition.insert(std::make_pair(vertexPositionsByIndex_[elements_[i][0]],     // order by position of left vertex
                                             std::make_pair(elements_[i], i)      // the element and its position in the insertion sequence
                                             ));

  auto it = std::get<0>(grid_->entityImps_[0]).begin();
  auto eIt = elementsByPosition.begin();

  // Looping over the vertices to get all elements assumes that the grid is connected
  for (size_t i=0; i<vertexPositions_.size()-1; i++, ++eIt)
  {
    OneDEntityImp<1,dimw,ct> newElement(0, grid_->getNextFreeId(), grid_->reversedBoundarySegmentNumbering_);
    newElement.vertex_[0] = it;
    it = it->succ_;
    newElement.vertex_[1] = it;

    newElement.levelIndex_ = eIt->second.second;
    newElement.leafIndex_  = eIt->second.second;

    grid_->elements(0).push_back(newElement);
  }

  // Create the index sets
  grid_->levelIndexSets_.resize(1);
  grid_->levelIndexSets_[0] = new OneDGridLevelIndexSet<const Grid>(*grid_, 0);
  grid_->levelIndexSets_[0]->setSizesAndTypes(vertexPositions_.size(), elements_.size());

  grid_->leafIndexSet_.setSizesAndTypes(vertexPositions_.size(), elements_.size());

  // Hand over the grid and delete the member pointer
  auto tmp = grid_;
  grid_ = nullptr;
  return std::unique_ptr<Grid>(tmp);
}

template <int dimw, class ct>
void GridFactory<OneDEmbeddedGrid<dimw,ct>>::
createBegin()
{
  vertexPositions_.clear();
}

template class GridFactory<OneDEmbeddedGrid<1,double>>;
template class GridFactory<OneDEmbeddedGrid<2,double>>;
template class GridFactory<OneDEmbeddedGrid<3,double>>;

} // namespace Dune
