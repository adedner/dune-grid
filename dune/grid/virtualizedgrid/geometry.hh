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

namespace Dune {

  template<int mydim, int coorddim, class GridImp>
  class VirtualizedGridGeometry :
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, VirtualizedGridGeometry>
  {
  public:
    typedef typename GridImp::ctype ctype;
    typedef ctype Volume;
    typedef FieldMatrix< ctype, mydim, coorddim > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddim, mydim > JacobianInverseTransposed;
    typedef FieldVector< ctype, coorddim > GlobalCoordinate;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual GeometryType type () const = 0;
      virtual bool affine() const = 0;
      virtual int corners () const = 0;
      virtual FieldVector<ctype, coorddim> corner (int i) const = 0;
      virtual FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const = 0;
      virtual JacobianTransposed jacobianTransposed ( const FieldVector<ctype, mydim>& local ) const = 0;
      virtual FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const = 0;
      virtual ctype integrationElement (const FieldVector<ctype, mydim>& local) const = 0;
      virtual JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const = 0;
    };

    template< class I >
    struct DUNE_PRIVATE Implementation final
      : public Interface
    {
      Implementation ( I&& i ) : impl_( std::forward<I>(i) ) {}
      Implementation *clone() const override { return new Implementation( *this ); }

      GeometryType type () const override { return impl().type(); }
      bool affine() const override { return impl().affine(); }
      int corners () const override { return impl().corners(); }
      FieldVector<ctype, coorddim> corner (int i) const override { return impl().corner(i); }
      FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const override { return impl().global(local); }
      JacobianTransposed jacobianTransposed ( const FieldVector<ctype, mydim>& local ) const override { return impl().jacobianTransposed(local); }
      FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const override { return impl().local(global); }
      ctype integrationElement (const FieldVector<ctype, mydim>& local) const override { return impl().integrationElement(local); }
      JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const override { return impl().jacobianInverseTransposed(local); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      I impl_;
    };
    // VIRTUALIZATION END

    struct Cache
    {
      template<class Geometry>
      void update (const Geometry& geometry)
      {
        type_ = geometry.type();
        affine_ = geometry.affine();
        corners_ = geometry.corners();

        corner_.resize(corners_);
        for (int i = 0; i < corners_; ++i)
          corner_[i] = geometry.corner(i);

        volume_ = geometry.volume();
        center_ = geometry.center();
      }

      //! return the element type identifier
      GeometryType type () const {
        return type_;
      }

      //! return wether we have an affine mapping
      bool affine () const {
        return affine_;
      }

      //! return the number of corners of this element. Corners are numbered 0...n-1
      int corners () const {
        return corners_;
      }

      //! access to coordinates of corners. Index is the number of the corner
      const FieldVector<ctype, coorddim>& corner (int i) const {
        return corner_[i];
      }

      //! return volume of the geometry
      Volume volume () const {
        return volume_;
      }

      //! return center of the geometry
      GlobalCoordinate center () const {
        return center_;
      }

    private:
      GeometryType type_{};
      bool affine_{};
      int corners_{};
      std::vector<FieldVector<ctype, coorddim>> corner_{};
      Volume volume_{};
      GlobalCoordinate center_{};
    };

  public:

    /** constructor from host geometry
     */
    template< class ImplGridGeometry >
    VirtualizedGridGeometry(ImplGridGeometry&& implGridGeometry)
      : impl_( new Implementation<ImplGridGeometry>( std::forward<ImplGridGeometry>(implGridGeometry) ) )
    {
      update(implGridGeometry);
    }

    VirtualizedGridGeometry(const VirtualizedGridGeometry& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    , cache_( other.cache_ )
    {}

    VirtualizedGridGeometry (VirtualizedGridGeometry&&) = default;

    VirtualizedGridGeometry& operator=(const VirtualizedGridGeometry& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
      cache_ = other.cache_;
      return *this;
    }

    VirtualizedGridGeometry& operator= (VirtualizedGridGeometry&&) = default;

    /** \brief Return the element type identifier
     */
    GeometryType type () const {
      return cache().type();
    }

    // return wether we have an affine mapping
    bool affine () const {
      return cache().affine();
    }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {
      return cache().corners();
    }

#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
    //! return volume of the geometry
    Volume volume () const {
      return cache().volume();
    }

    //! return center of the geometry
    GlobalCoordinate center () const {
      return cache().center();
    }
#endif

    //! access to coordinates of corners. Index is the number of the corner
    FieldVector<ctype, coorddim> corner (int i) const {
      return cache().corner(i);
    }

    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const {
      return impl_->global(local);
    }

    /** \brief Return the transposed of the Jacobian
     */
    JacobianTransposed
    jacobianTransposed ( const FieldVector<ctype, mydim>& local ) const {
      return impl_->jacobianTransposed(local);
    }

    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const {
      return impl_->local(global);
    }

    ctype integrationElement (const FieldVector<ctype, mydim>& local) const {
      return impl_->integrationElement(local);
    }
    //! The Jacobian matrix of the mapping from the reference element to this element
    JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
      return impl_->jacobianInverseTransposed(local);
    }

    const auto& cache() const {
#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
      return cache_;
#else
      return *impl_;
#endif
    }

    template<class Grid>
    void update (const Grid& grid) {
#ifndef DUNE_VIRTUALIZEDGRID_NO_CACHE
      cache_.update(grid);
#endif
    }

    std::unique_ptr<Interface> impl_;
    Cache cache_;
  };

}  // namespace Dune

#endif
