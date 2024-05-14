// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <dune/common/reservedvector.hh>
#include <dune/grid/onedgrid/onedgridfactory.hh>
#include <dune/grid/onedgrid/onedgridindexsets.hh>

namespace Dune {

template <int dimw, class ct>
GridFactory<OneDEmbeddedGrid<dimw,ct>>::
GridFactory() :
  factoryOwnsGrid_(true)
{
  grid_ = new OneDEmbeddedGrid<dimw,ct>;

  createBegin();
}

template <int dimw, class ct>
GridFactory<OneDEmbeddedGrid<dimw,ct>>::
GridFactory(Grid* grid) :
  factoryOwnsGrid_(false)
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
  vertexPositions_.push_back(pos);
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
  std::size_t idx = intersection.inside().impl().leafIndex();
  for(const auto& bidx : boundarySegments_)
    if(idx == bidx)
      return true;
  return false;
}

template <int dimw, class ct>
unsigned int GridFactory<OneDEmbeddedGrid<dimw,ct>>::
insertionIndex(const typename Grid::LeafIntersection& intersection) const
{
  std::size_t idx = intersection.inside().impl().leafIndex();
  for(std::size_t i = 0; i < boundarySegments_.size(); ++i)
    if(idx == boundarySegments_[i])
      return i;
  return (unsigned int)(-1);
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

  std::vector<Dune::ReservedVector<std::size_t,2>> vertexToElements(vertexPositions_.size());
  for (std::size_t i = 0; i < elements_.size(); ++i) {
    auto [v0,v1] = elements_[i];
    vertexToElements[v0].push_back(i);
    vertexToElements[v1].push_back(i);
  }

  // Check the numbering of the boundary segments
  if (!(boundarySegments_.empty() || boundarySegments_.size() == 2))
    DUNE_THROW(GridError, "You cannot provide more than two boundary segments to a OneDGrid (it must be connected).");

  // order the elements, such that the 2nd vertex of the ith element corresponds to
  // the 1st vertex of the (i+1)th element.
  std::vector<std::pair<std::array<unsigned int, 2>,std::size_t>> sortedElements;
  sortedElements.reserve(elements_.size());

  // The first vertex is either the first boundarySegment index or first first element vertex
  std::size_t v0 = boundarySegments_.size() > 0 ? boundarySegments_[0] : elements_[0][0];
  std::size_t v1 = v0;
  std::size_t e1 = std::size_t(-1);
  while (sortedElements.size() < elements_.size()) {
    for (std::size_t e : vertexToElements[v1]) {
      if (e != e1) {
        if (elements_[e][0] != v1)
          std::swap(elements_[e][0], elements_[e][1]);

        sortedElements.push_back(std::make_pair(elements_[e],e));
        v1 = elements_[e][1];
        e1 = e;
        break;
      }
    }
  }

  // Insert the vertices into the grid
  grid_->entityImps_.resize(1);

  v0 = sortedElements[0].first[0];

  { // insert the first vertex into the grid
    OneDEntityImp<0,dimw,ct> newVertex(0, vertexPositions_[v0], grid_->getNextFreeId());
    newVertex.leafIndex_ = v0;
    newVertex.levelIndex_ = v0;
    grid_->vertices(0).push_back(newVertex);
  }

  // insert the second vertex of all elements in order of the elements
  for (const auto& e : sortedElements)
  {
    auto vIt = grid_->vertices(0).rbegin();

    v1 = e.first[1];
    if (v1 != v0) {
      OneDEntityImp<0,dimw,ct> newVertex(0, vertexPositions_[v1], grid_->getNextFreeId());
      newVertex.leafIndex_ = v1;
      newVertex.levelIndex_ = v1;
      grid_->vertices(0).push_back(newVertex);
    } else {
      // connect the first and the last vertex
      auto rbegin = grid_->vertices(0).rbegin();
      auto begin = grid_->vertices(0).begin();
      rbegin->succ_ = begin;
      begin->pred_ = rbegin;
    }

    OneDEntityImp<1,dimw,ct> newElement(0, grid_->getNextFreeId(), false);
    newElement.vertex_[0] = vIt;
    vIt = vIt->succ_;
    newElement.vertex_[1] = vIt;
    newElement.levelIndex_ = e.second;
    newElement.leafIndex_  = e.second;
    grid_->elements(0).push_back(newElement);

    // periodic connectivity of the elements
    if (v0 == v1) {
      auto rbegin = grid_->elements(0).rbegin();
      auto begin = grid_->elements(0).begin();
      rbegin->succ_ = begin;
      begin->pred_ = rbegin;
    }
  }

  if (v0 == v1 && boundarySegments_.size() > 0)
    DUNE_THROW(GridError, "Boundary segments given for a periodic grid!");

  if (boundarySegments_.size() > 1 && (v0 != boundarySegments_[0] || v1 != boundarySegments_[1])) {
    DUNE_THROW(GridError, "The boundary segments are not at the boundary!");
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
