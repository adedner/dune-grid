// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONE_D_GRID_INTERSECTIONIT_HH
#define DUNE_ONE_D_GRID_INTERSECTIONIT_HH

/** \file
 * \brief The OneDGridIntersectionIterator class
 */

namespace Dune {

  //**********************************************************************
  //
  // --OneDGridIntersectionIterator
  // --IntersectionIterator
  /** \brief Iterator over all element neighbors
   * \ingroup OneDGrid
     Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
     a neighbor is an entity of codimension 0 which has a common entity of codimension
     These neighbors are accessed via a IntersectionIterator. This allows the implement
     non-matching meshes. The number of neigbors may be different from the number
     of an element!
   */
  template<class GridImp>
  class OneDGridLevelIntersectionIterator
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };

    friend class OneDGridEntity<0,dim,GridImp>;

    //! Constructor for a given grid entity and a given neighbor
    OneDGridLevelIntersectionIterator(OneDEntityImp<1>* center, int nb)
      : center_(center), neighbor_(nb),
        intersectionSelfLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionNeighborLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionGlobal_(OneDGridGeometry<0,1,GridImp>())
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    OneDGridLevelIntersectionIterator(OneDEntityImp<1>* center)
      : center_(center), neighbor_(2),
        intersectionSelfLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionNeighborLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionGlobal_(OneDGridGeometry<0,1,GridImp>())
    {}

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef Dune::Intersection<const GridImp, Dune::OneDGridLevelIntersectionIterator> Intersection;

    //! equality
    bool equals(const OneDGridLevelIntersectionIterator<GridImp>& other) const {
      return (center_ == other.center_) && (neighbor_ == other.neighbor_);
    }

    //! prefix increment
    void increment() {
      neighbor_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      return reinterpret_cast<const Intersection&>(*this);
    }

    OneDEntityImp<1>* target() const {
      const bool isValid = center_ && neighbor_>=0 && neighbor_<2;

      if (!isValid)
        return center_;
      else if (neighbor_==0)
        return center_->pred_;
      else
        return center_->succ_;

    }

    //! return true if intersection is with boundary.
    bool boundary () const {

      // Check whether we're on the left boundary
      if (neighbor_==0) {

        // If there's an element to the left we can't be on the boundary
        if (center_->pred_)
          return false;

        const OneDEntityImp<1>* ancestor = center_;

        while (ancestor->level_!=0) {

          // Check if we're the left son of our father
          if (ancestor != ancestor->father_->sons_[0])
            return false;

          ancestor = ancestor->father_;
        }

        // We have reached level 0.  If there is no element of the left
        // we're truly on the boundary
        return !ancestor->pred_;
      }

      // ////////////////////////////////
      //   Same for the right boundary
      // ////////////////////////////////
      // If there's an element to the right we can't be on the boundary
      if (center_->succ_)
        return false;

      const OneDEntityImp<1>* ancestor = center_;

      while (ancestor->level_!=0) {

        // Check if we're the left son of our father
        if (ancestor != ancestor->father_->sons_[1])
          return false;

        ancestor = ancestor->father_;
      }

      // We have reached level 0.  If there is no element of the left
      // we're truly on the boundary
      return !ancestor->succ_;

    }

    //! return true if across the edge a neighbor on this level exists
    bool neighbor () const {
      assert(neighbor_ >= 0 && neighbor_ < 2);

      return (neighbor_==0)
             ? center_->pred_ && center_->pred_->vertex_[1] == center_->vertex_[0]
             : center_->succ_ && center_->succ_->vertex_[0] == center_->vertex_[1];

    }

    //! return true if intersection is conform.
    bool conforming () const {
      return true;
    }

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const
    {
      return OneDGridEntityPointer<0,GridImp>(center_);
    }

    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const
    {
      assert(neighbor());
      return OneDGridEntityPointer<0,GridImp>(target());
    }

    //! return information about the Boundary
    int boundaryId () const {
      return 1;
    }

    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &geometryInInside () const
    {
      GridImp::getRealImplementation(intersectionSelfLocal_).setPosition( (indexInInside() == 0) ? 0 : 1 );
      return intersectionSelfLocal_;
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry &geometryInOutside () const
    {
      GridImp::getRealImplementation(intersectionNeighborLocal_).setPosition( (indexInInside() == 0) ? 1 : 0 );
      return intersectionNeighborLocal_;
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    const Geometry &geometry () const
    {
      GridImp::getRealImplementation(intersectionGlobal_).target_ = center_->vertex_[neighbor_];
      return intersectionGlobal_;
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return GeometryType( 0 );
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return neighbor_;
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const
    {
      // If numberInSelf is 0 then numberInNeighbor is 1 and vice versa
      return 1-neighbor_;
    }

    //! return outer normal
    const FieldVector<typename GridImp::ctype, dimworld>& outerNormal (const FieldVector<typename GridImp::ctype, dim-1>& local) const {
      outerNormal_[0] = (neighbor_==0) ? -1 : 1;
      return outerNormal_;
    }

    //! Return outer normal scaled with the integration element
    const FieldVector<typename GridImp::ctype, dimworld>& integrationOuterNormal (const FieldVector<typename GridImp::ctype, dim-1>& local) const
    {
      return outerNormal(local);
    }

    //! return unit outer normal
    const FieldVector<typename GridImp::ctype, dimworld>& unitOuterNormal (const FieldVector<typename GridImp::ctype, dim-1>& local) const {
      return outerNormal(local);
    }

  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    OneDEntityImp<1>* center_;

    //! vector storing the outer normal
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;

    /** \brief Count on which neighbor we are lookin' at.  Can be only 0 or 1. */
    int neighbor_;

    /** \brief The geometry that's being returned when intersectionSelfLocal() is called
     */
    mutable MakeableInterfaceObject<LocalGeometry> intersectionSelfLocal_;

    /** \brief The geometry that's being returned when intersectionNeighborLocal() is called
     */
    mutable MakeableInterfaceObject<LocalGeometry> intersectionNeighborLocal_;

    //! The geometry that's being returned when intersectionSelfGlobal() is called
    mutable MakeableInterfaceObject<Geometry> intersectionGlobal_;

  };


  template<class GridImp>
  class OneDGridLeafIntersectionIterator
  {
    enum { dim=GridImp::dimension };
    enum { dimworld=GridImp::dimensionworld };

    friend class OneDGridEntity<0,dim,GridImp>;

    //! Constructor for a given grid entity and a given neighbor
    OneDGridLeafIntersectionIterator(OneDEntityImp<1>* center, int nb)
      : center_(center), neighbor_(nb),
        intersectionSelfLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionNeighborLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionGlobal_(OneDGridGeometry<0,1,GridImp>())
    {}

    /** \brief Constructor creating the 'one-after-last'-iterator */
    OneDGridLeafIntersectionIterator(OneDEntityImp<1>* center)
      : center_(center), neighbor_(2),
        intersectionSelfLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionNeighborLocal_(OneDGridGeometry<0,1,GridImp>()),
        intersectionGlobal_(OneDGridGeometry<0,1,GridImp>())
    {}

  public:

    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef Dune::Intersection<const GridImp, Dune::OneDGridLeafIntersectionIterator> Intersection;

    //! The Destructor
    ~OneDGridLeafIntersectionIterator() {};

    //! equality
    bool equals(const OneDGridLeafIntersectionIterator<GridImp>& other) const {
      return (center_ == other.center_) && (neighbor_ == other.neighbor_);
    }

    //! prefix increment
    void increment() {
      neighbor_++;
    }

    //! \brief dereferencing
    const Intersection & dereference() const
    {
      return reinterpret_cast<const Intersection&>(*this);
    }

    OneDEntityImp<1>* target() const {
      const bool isValid = center_ && neighbor_>=0 && neighbor_<2;

      if (!isValid)
        return center_;

      if (neighbor_==0) {

        // Get left leaf neighbor
        if (center_->pred_ && center_->pred_->vertex_[1] == center_->vertex_[0]) {

          OneDEntityImp<1>* leftLeafNeighbor = center_->pred_;
          while (!leftLeafNeighbor->isLeaf()) {
            assert (leftLeafNeighbor->sons_[1] != NULL);
            leftLeafNeighbor = leftLeafNeighbor->sons_[1];
          }
          return leftLeafNeighbor;

        } else {

          OneDEntityImp<1>* ancestor = center_;
          while (ancestor->father_) {
            ancestor = ancestor->father_;
            if (ancestor->pred_ && ancestor->pred_->vertex_[1] == ancestor->vertex_[0]) {
              assert(ancestor->pred_->isLeaf());
              return ancestor->pred_;
            }
          }

          DUNE_THROW(GridError, "Programming error, apparently we're on the left boundary, neighbor_==2 should not occur!");
        }

      } else {

        // Get right leaf neighbor
        if (center_->succ_ && center_->succ_->vertex_[0] == center_->vertex_[1]) {

          OneDEntityImp<1>* rightLeafNeighbor = center_->succ_;
          while (!rightLeafNeighbor->isLeaf()) {
            assert (rightLeafNeighbor->sons_[0] != NULL);
            rightLeafNeighbor = rightLeafNeighbor->sons_[0];
          }
          return rightLeafNeighbor;

        } else {

          OneDEntityImp<1>* ancestor = center_;
          while (ancestor->father_) {
            ancestor = ancestor->father_;
            if (ancestor->succ_ && ancestor->succ_->vertex_[0] == ancestor->vertex_[1]) {
              assert(ancestor->succ_->isLeaf());
              return ancestor->succ_;
            }
          }

          DUNE_THROW(GridError, "Programming error, apparently we're on the right boundary, neighbor_==3 should not occur!");
        }

      }

    }

    //! return true if intersection is with boundary.
    bool boundary () const {

      // Check whether we're on the left boundary
      if (neighbor_==0) {

        // If there's an element to the left we can't be on the boundary
        if (center_->pred_)
          return false;

        const OneDEntityImp<1>* ancestor = center_;

        while (ancestor->level_!=0) {

          // Check if we're the left son of our father
          if (ancestor != ancestor->father_->sons_[0])
            return false;

          ancestor = ancestor->father_;
        }

        // We have reached level 0.  If there is no element of the left
        // we're truly on the boundary
        return !ancestor->pred_;
      }

      // ////////////////////////////////
      //   Same for the right boundary
      // ////////////////////////////////
      // If there's an element to the right we can't be on the boundary
      if (center_->succ_)
        return false;

      const OneDEntityImp<1>* ancestor = center_;

      while (ancestor->level_!=0) {

        // Check if we're the left son of our father
        if (ancestor != ancestor->father_->sons_[1])
          return false;

        ancestor = ancestor->father_;
      }

      // We have reached level 0.  If there is no element of the left
      // we're truly on the boundary
      return !ancestor->succ_;

    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor () const {
      return !boundary();
    }

    //! return true if intersection is conform.
    bool conforming () const {
      return true;
    }

    //! return EntityPointer to the Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    EntityPointer inside() const
    {
      return OneDGridEntityPointer<0,GridImp>(center_);
    }

    //! return EntityPointer to the Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    EntityPointer outside() const
    {
      return OneDGridEntityPointer<0,GridImp>(target());
    }

    //! return information about the Boundary
    int boundaryId () const {
      return 1;
    }

    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    const LocalGeometry &geometryInInside () const
    {
      GridImp::getRealImplementation(intersectionSelfLocal_).setPosition( (indexInInside() == 0) ? 0 : 1 );
      return intersectionSelfLocal_;
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    const LocalGeometry &geometryInOutside () const
    {
      GridImp::getRealImplementation(intersectionNeighborLocal_).setPosition( (indexInInside() == 0) ? 1 : 0 );
      return intersectionNeighborLocal_;
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    const Geometry &geometry () const
    {
      GridImp::getRealImplementation(intersectionGlobal_).target_ = center_->vertex_[neighbor_%2];
      return intersectionGlobal_;
    }

    /** \brief obtain the type of reference element for this intersection */
    GeometryType type () const
    {
      return GeometryType( 0 );
    }

    //! local index of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return neighbor_ % 2;
    }

    //! local index of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const
    {
      // If numberInSelf is 0 then numberInNeighbor is 1 and vice versa
      return 1-(neighbor_ % 2);
    }

    //! return outer normal
    const FieldVector<typename GridImp::ctype, dimworld>& outerNormal (const FieldVector<typename GridImp::ctype, dim-1>& local) const {
      outerNormal_[0] = ((neighbor_%2)==0) ? -1 : 1;
      return outerNormal_;
    }

    //! Return outer normal scaled with the integration element
    const FieldVector<typename GridImp::ctype, dimworld>& integrationOuterNormal (const FieldVector<typename GridImp::ctype, dim-1>& local) const
    {
      return outerNormal(local);
    }

    //! return unit outer normal
    const FieldVector<typename GridImp::ctype, dimworld>& unitOuterNormal (const FieldVector<typename GridImp::ctype, dim-1>& local) const {
      return outerNormal(local);
    }

  private:
    //**********************************************************
    //  private methods
    //**********************************************************

    OneDEntityImp<1>* center_;

    //! vector storing the outer normal
    mutable FieldVector<typename GridImp::ctype, dimworld> outerNormal_;

    /** \brief Count on which neighbor we are lookin' at

       0,1 are the level neighbors, 2 and 3 are the leaf neighbors,
       if they differ from the level neighbors. */
    int neighbor_;

    /** \brief The geometry that's being returned when intersectionSelfLocal() is called
     */
    mutable MakeableInterfaceObject<LocalGeometry> intersectionSelfLocal_;

    /** \brief The geometry that's being returned when intersectionNeighborLocal() is called
     */
    mutable MakeableInterfaceObject<LocalGeometry> intersectionNeighborLocal_;

    //! The geometry that's being returned when intersectionSelfGlobal() is called
    mutable MakeableInterfaceObject<Geometry> intersectionGlobal_;

  };

}  // namespace Dune

#endif
