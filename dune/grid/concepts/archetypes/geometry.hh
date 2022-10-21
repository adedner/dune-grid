// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_CONCEPTS_ARCHETYPES_GEOMETRY_HH
#define DUNE_GRID_CONCEPTS_ARCHETYPES_GEOMETRY_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

namespace Dune::Concept::Archetypes {

struct ReferenceElement {};

template <int mydim, int cdim = mydim>
struct Geometry
{
  static constexpr int mydimension = mydim;
  static constexpr int coorddimension = cdim;

  using ctype = double;
  using LocalCoordinate = Dune::FieldVector<ctype, mydim>;
  using GlobalCoordinate = Dune::FieldVector<ctype, cdim>;
  using Volume = ctype;
  using Jacobian = Dune::FieldMatrix<ctype, cdim, mydim>;
  using JacobianTransposed = Dune::FieldMatrix<ctype, mydim, cdim>;
  using JacobianInverse = Dune::FieldMatrix<ctype, mydim, cdim>;
  using JacobianInverseTransposed = Dune::FieldMatrix<ctype, cdim, mydim>;

  Dune::GeometryType type () const;
  bool affine () const;
  int corners() const;

  GlobalCoordinate corner(int i) const;
  GlobalCoordinate global(const LocalCoordinate&) const;
  LocalCoordinate local(const GlobalCoordinate&) const;
  GlobalCoordinate center() const;

  Volume integrationElement(const LocalCoordinate&) const;
  Volume volume() const;

  Jacobian jacobian(const LocalCoordinate&) const;
  JacobianTransposed jacobianTransposed(const LocalCoordinate&) const;
  JacobianInverse jacobianInverse(const LocalCoordinate&) const;
  JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate&) const;
};

template <int mydim, int cdim>
Archetypes::ReferenceElement referenceElement(const Geometry<mydim,cdim>&);

} // end namespace Dune::Concept::Archetypes


#endif // DUNE_GRID_CONCEPTS_ARCHETYPES_GEOMETRY_HH
