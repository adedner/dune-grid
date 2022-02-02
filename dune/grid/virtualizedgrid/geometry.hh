// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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
    typedef FieldMatrix< ctype, mydim, coorddim > JacobianTransposed;
    typedef FieldMatrix< ctype, coorddim, mydim > JacobianInverseTransposed;

  private:
    // VIRTUALIZATION BEGIN
    struct Interface
    {
      virtual ~Interface () = default;
      virtual Interface *clone () const = 0;
      virtual GeometryType type () const = 0;
      virtual bool affine() const = 0;
      virtual int corners () const = 0;
      virtual const FieldVector<ctype, coorddim> corner (int i) const = 0;
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
      Implementation ( const I& i ) : impl_( i ) {}
      virtual Implementation *clone() const override { return new Implementation( *this ); }

      virtual GeometryType type () const override { return impl().type(); }
      virtual bool affine() const override { return impl().affine(); }
      virtual int corners () const override { return impl().corners(); }
      virtual const FieldVector<ctype, coorddim> corner (int i) const override { return impl().corner(i); }
      virtual FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const override { return impl().global(local); }
      virtual JacobianTransposed jacobianTransposed ( const FieldVector<ctype, mydim>& local ) const override { return impl().jacobianTransposed(local); }
      virtual FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const override { return impl().local(global); }
      virtual ctype integrationElement (const FieldVector<ctype, mydim>& local) const override { return impl().integrationElement(local); }
      virtual JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const override { return impl().jacobianInverseTransposed(local); }

    private:
      const auto &impl () const { return impl_; }
      auto &impl () { return impl_; }

      const I& impl_;
    };
    // VIRTUALIZATION END

  public:

    /** constructor from host geometry
     */
    template< class ImplGridGeometry >
    VirtualizedGridGeometry(const ImplGridGeometry& implGridGeometry)
      : impl_( new Implementation<ImplGridGeometry>(implGridGeometry) )
    {}

    VirtualizedGridGeometry(const VirtualizedGridGeometry& other)
    : impl_( other.impl_ ? other.impl_->clone() : nullptr )
    {}

    VirtualizedGridGeometry ( VirtualizedGridGeometry && ) = default;

    VirtualizedGridGeometry& operator=(const VirtualizedGridGeometry& other)
    {
      impl_.reset( other.impl_ ? other.impl_->clone() : nullptr );
    }

    /** \brief Return the element type identifier
     */
    GeometryType type () const {
      return impl_->type();
    }

    // return wether we have an affine mapping
    bool affine() const {
      return impl_->affine();
    }

    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {
      return impl_->corners();
    }


    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<ctype, coorddim> corner (int i) const {
      return impl_->corner(i);
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


    /**
     */
    ctype integrationElement (const FieldVector<ctype, mydim>& local) const {
      return impl_->integrationElement(local);
    }


    //! The Jacobian matrix of the mapping from the reference element to this element
    JacobianInverseTransposed jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
      return impl_->jacobianInverseTransposed(local);
    }


    std::unique_ptr<Interface> impl_;

  };

}  // namespace Dune

#endif
