// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_VIRTUALIZEDGRIDGEOMETRY_HH
#define DUNE_VIRTUALIZEDGRIDGEOMETRY_HH

/** \file
 * \brief The VirtualizedGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/virtualizedgrid/common/typeerasure.hh>

namespace Dune {

  template<int mydim, int coorddim, class GridImp>
  struct VirtualizedGridGeometryDefinition
  {
    using ctype = typename GridImp::ctype;
    using Jacobian = FieldMatrix<ctype, coorddim, mydim>;
    using JacobianTransposed = FieldMatrix<ctype, mydim, coorddim>;
    using JacobianInverse = FieldMatrix<ctype, mydim, coorddim>;
    using JacobianInverseTransposed = FieldMatrix<ctype, coorddim, mydim>;

    struct Interface
    {
      virtual ~Interface () = default;
      virtual GeometryType type () const = 0;
      virtual bool affine () const = 0;
      virtual int corners () const = 0;
      virtual FieldVector<ctype, coorddim> corner (int i) const = 0;
      virtual FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const = 0;
      virtual FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const = 0;
      virtual ctype integrationElement (const FieldVector<ctype, mydim>& local) const = 0;
      virtual Jacobian jacobian (const FieldVector<ctype, mydim>& local) const = 0;
      virtual JacobianTransposed jacobianTransposed (const FieldVector<ctype, mydim>& local) const = 0;
      virtual JacobianInverse jacobianInverse (const FieldVector<ctype, mydim>& local) const = 0;
      virtual JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const = 0;
    };

    template<class Wrapper>
    struct Implementation
      : public Wrapper
    {
      using Wrapper::Wrapper;

      GeometryType type () const final { return this->get().type(); }
      bool affine () const final { return this->get().affine(); }
      int corners () const final { return this->get().corners(); }
      FieldVector<ctype, coorddim> corner (int i) const final { return this->get().corner(i); }
      FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const final { return this->get().global(local); }
      FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const final { return this->get().local(global); }
      ctype integrationElement (const FieldVector<ctype, mydim>& local) const final { return this->get().integrationElement(local); }
      Jacobian jacobian (const FieldVector<ctype, mydim>& local) const final { return this->get().jacobian(local); }
      JacobianTransposed jacobianTransposed (const FieldVector<ctype, mydim>& local) const final { return this->get().jacobianTransposed(local); }
      JacobianInverse jacobianInverse (const FieldVector<ctype, mydim>& local) const final { return this->get().jacobianInverse(local); }
      JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const final { return this->get().jacobianInverseTransposed(local); }
    };

    using Base = Polymorphic::TypeErasureBase<Interface, Implementation>;
  };


  template<int mydim, int coorddim, class GridImp>
  class VirtualizedGridGeometry :
    public VirtualizedGridGeometryDefinition<mydim,coorddim,GridImp>::Base,
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, VirtualizedGridGeometry>
  {
    using Definition = VirtualizedGridGeometryDefinition<mydim,coorddim,GridImp>;
    using Base = typename Definition::Base;

  public:
    using ctype = typename GridImp::ctype;
    using Jacobian = FieldMatrix<ctype, coorddim, mydim>;
    using JacobianTransposed = FieldMatrix<ctype, mydim, coorddim>;
    using JacobianInverse = FieldMatrix<ctype, mydim, coorddim>;
    using JacobianInverseTransposed = FieldMatrix<ctype, coorddim, mydim>;

  public:

    /** constructor from host geometry
     */
    template <class Impl, disableCopyMove<VirtualizedGridGeometry,Impl> = 0>
    VirtualizedGridGeometry (Impl&& impl)
      : Base{std::forward<Impl>(impl)}
    {}

    /** \brief Return the element type identifier
     */
    GeometryType type () const
    {
      return this->asInterface().type();
    }

    // return wether we have an affine mapping
    bool affine () const
    {
      return this->asInterface().affine();
    }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const
    {
      return this->asInterface().corners();
    }

    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<ctype, coorddim> corner (int i) const
    {
      return this->asInterface().corner(i);
    }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const
    {
      return this->asInterface().global(local);
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const
    {
      return this->asInterface().local(global);
    }

    ctype integrationElement (const FieldVector<ctype, mydim>& local) const
    {
      return this->asInterface().integrationElement(local);
    }

    /** \brief Return the transposed of the Jacobian
     */
    Jacobian jacobian (const FieldVector<ctype, mydim>& local) const
    {
      return this->asInterface().jacobian(local);
    }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed jacobianTransposed (const FieldVector<ctype, mydim>& local) const
    {
      return this->asInterface().jacobianTransposed(local);
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    JacobianInverse jacobianInverse (const FieldVector<ctype, mydim>& local) const
    {
      return this->asInterface().jacobianInverse(local);
    }

    //! The Jacobian matrix of the mapping from the reference element to this element
    JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const
    {
      return this->asInterface().jacobianInverseTransposed(local);
    }
  };

}  // namespace Dune

#endif
