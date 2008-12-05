// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \file
 * \brief Encapsulates some UG macros and functions
 */

/** \todo Here only to provide the constant DBL_EPSILON.  There's maybe a better way? */
#include "float.h"

#if defined ModelP && UG_DIM == 3

/* This is a hack.  The following define is from the UG header pargm.h.  We can include that
   header only once, because it contains dimension-independent stuff.  However, when we get
   here for the second time for the 3d stuff, the macro PRIO2LISTPART has been undefined.
   Hence we define it again here

   /* define mapping from object priority to position in linked list */
#define PRIO2LISTPART(listtype,prio)                                         \
  ((listtype == ELEMENT_LIST) ? ((prio == PrioHGhost) ? 0 :                \
                                 (prio == PrioVGhost) ? 0 : (prio == PrioVHGhost) ? 0 :               \
                                 (prio == PrioMaster) ? 1 : -1) :                                 \
   ((prio == PrioHGhost) ? 0 : (prio ==PrioVGhost) ? 0 :            \
      (prio == PrioVHGhost) ? 0 :                                  \
      (prio == PrioBorder) ? 2 : (prio == PrioMaster) ? 2 : -1))
#endif

namespace Dune {

  /** \brief Encapsulates a few UG methods and macros
   *
   * This class provides a wrapper to several methods and macros from
   * UG.  There are two reasons for doing this.  First, we don't want
   * to call UG macros directly from DUNE, because they pollute the
   * namespace and therefore we undefine them all.  Secondly,  UG methods
   * appear in the namespaces UG::D2 and UG::D3, but we need the dimension
   * as a template parameter.
   */
#if UG_DIM == 2
  template<int dim>
  class UG_NS {};
#endif

  template<>
  class UG_NS< UG_DIM > {
  public:

    // //////////////////////////////////////////////
    //   Types exported by UG
    // //////////////////////////////////////////////

#if UG_DIM == 2
#define UG_NAMESPACE UG::D2
#else
#define UG_NAMESPACE UG::D3
#endif

    typedef UG_NAMESPACE ::RefinementRule RefinementRule;

    typedef UG_NAMESPACE ::CoeffProcPtr CoeffProcPtr;

    typedef UG_NAMESPACE ::UserProcPtr UserProcPtr;

#ifndef UG_LGMDOMAIN
    typedef UG_NAMESPACE ::BndSegFuncPtr BndSegFuncPtr;
#endif

    /** \todo This type becomes obsolete as soon as UG implements faces and edges */
    typedef UG_NAMESPACE ::vector Vector;

    typedef UG_NAMESPACE ::multigrid MultiGrid;

    typedef UG_NAMESPACE ::grid Grid;

    typedef UG_NAMESPACE ::edge Edge;

    typedef UG_NAMESPACE ::node Node;

    typedef UG_NAMESPACE ::element Element;

    typedef UG_NAMESPACE ::vertex Vertex;

    /** \brief Types of the subentities parametrized by the codimension.  Gets specialized below */
    template <int codim>
    class Entity;

    // //////////////////////////////////////////////
    //   Constants exported by UG
    // //////////////////////////////////////////////

    enum {GM_REFINE_NOT_CLOSED = UG_NAMESPACE ::GM_REFINE_NOT_CLOSED};

    enum {GM_COPY_ALL = UG_NAMESPACE ::GM_COPY_ALL};

    enum {GM_REFINE_TRULY_LOCAL = UG_NAMESPACE ::GM_REFINE_TRULY_LOCAL};

    enum {GM_REFINE_PARALLEL = UG_NAMESPACE ::GM_REFINE_PARALLEL};

    enum {GM_REFINE_NOHEAPTEST = UG_NAMESPACE ::GM_REFINE_NOHEAPTEST};

    /** \brief Control word entries */
    enum {NEWEL_CE      = UG_NAMESPACE ::NEWEL_CE,
          COARSEN_CE    = UG_NAMESPACE ::COARSEN_CE,
          ECLASS_CE     = UG_NAMESPACE ::ECLASS_CE,
          MARK_CE       = UG_NAMESPACE ::MARK_CE};

    /** \brief Refinement rules */
    enum {NO_REFINEMENT = UG_NAMESPACE ::NO_REFINEMENT,
          RED           = UG_NAMESPACE ::RED,
          COARSE        = UG_NAMESPACE ::COARSE};

    enum {RED_CLASS = UG_NAMESPACE ::RED_CLASS};

    enum {GM_OK = UG_NAMESPACE ::GM_OK};

    enum {MAX_SONS = UG_NAMESPACE ::MAX_SONS};

    /** \brief The PFIRSTNODE macro which returns the first node in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Node* PFirstNode(const UG_NS< UG_DIM >::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::ELEMENT_LIST;
      using UG::NODE_LIST;
      return PFIRSTNODE(grid);
    }

    /** \brief The FIRSTNODE macro which returns the first node in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Node* FirstNode(UG_NS< UG_DIM >::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::ELEMENT_LIST;
      using UG::NODE_LIST;
      return FIRSTNODE(grid);
    }

    /** \brief The PFIRSTELEMENT macro which returns the first element in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Element* PFirstElement(const UG_NS< UG_DIM >::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::ELEMENT_LIST;
      using UG::NODE_LIST;
      return PFIRSTELEMENT(grid);
    }

    /** \brief The FIRSTELEMENT macro which returns the first element in a
     * grid even in a parallel setting.
     */
    static UG_NS< UG_DIM >::Element* FirstElement(UG_NS< UG_DIM >::Grid* grid) {
      using UG::PrioHGhost;
      using UG::PrioVGhost;
      using UG::PrioVHGhost;
      using UG::PrioMaster;
      using UG::PrioBorder;
      using UG::ELEMENT_LIST;
      return FIRSTELEMENT(grid);
    }

    /** \brief Returns pointers to the coordinate arrays of an UG element */
    static void Corner_Coordinates(const UG_NS< UG_DIM >::Element* theElement, double* x[]) {
      using UG_NAMESPACE ::NODE;
      using UG_NAMESPACE ::TRIANGLE;
      using UG_NAMESPACE ::QUADRILATERAL;
      using UG_NAMESPACE ::TETRAHEDRON;
      using UG_NAMESPACE ::PYRAMID;
      using UG_NAMESPACE ::PRISM;
      using UG_NAMESPACE ::n_offset;
      using UG::UINT;
      int n;    // Dummy variable just to please the macro
      CORNER_COORDINATES(theElement, n, x);
    }

    static int GlobalToLocal(int n, const double** cornerCoords,
                             const double* EvalPoint, double* localCoord) {
      return UG_NAMESPACE ::UG_GlobalToLocal(n, cornerCoords, EvalPoint, localCoord);
    }

    /** \brief Computes the element volume */
    static double Area_Of_Element(int n, const double** cornerCoords) {
      double area;
      using UG::DOUBLE;
      using UG_NAMESPACE ::DOUBLE_VECTOR;
#if UG_DIM == 2
      AREA_OF_ELEMENT_2D(n,cornerCoords,area);
#else
      AREA_OF_ELEMENT_3D(n,cornerCoords,area);
#endif
      return area;
    }

    static int myLevel (const UG_NS< UG_DIM >::Element* theElement) {
      using UG::UINT;
      return LEVEL(theElement);
    }

    static int myLevel (const UG_NS< UG_DIM >::Node* theNode) {
      using UG::UINT;
      return LEVEL(theNode);
    }

    //! return true if element has an exact copy on the next level
    static bool hasCopy (const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::control_entries;
      using UG::UINT;
      using UG_NAMESPACE ::REFINECLASS_CE;
      using UG_NAMESPACE ::YELLOW_CLASS;
      return REFINECLASS(theElement) == YELLOW_CLASS;
    }

    //! return true if element has an exact copy on the next level
    static bool isRegular (const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::control_entries;
      using UG::UINT;
      return ECLASS(theElement) == RED_CLASS;
    }

    //! \todo Please doc me!
    static int Sides_Of_Elem(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return SIDES_OF_ELEM(theElement);
    }

    //! Encapsulates the NBELEM macro
    static UG_NS<UG_DIM>::Element* NbElem(const UG_NS< UG_DIM >::Element* theElement, int nb) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::nb_offset;
      using UG::UINT;
      return NBELEM(theElement, nb);
    }

    static int boundaryId(const UG_NS< UG_DIM >::Element* theElement, int nb) {
      using UG_NAMESPACE ::BNDS;
      using UG::UINT;
      using UG_NAMESPACE ::side_offset;

      BNDS* bnds = ELEM_BNDS(theElement,nb);
      int id = UG_NAMESPACE ::GetBoundarySegmentId(bnds);
      return id;
    }

    //! Returns true if the i-th side of the element is on the domain boundary
    static bool Side_On_Bnd(const UG_NS< UG_DIM >::Element* theElement, int i) {
      using UG_NAMESPACE ::BNDS;
      using UG_NAMESPACE ::BEOBJ;
      using UG_NAMESPACE ::side_offset;
      using UG::UINT;
      using UG_NAMESPACE ::GM_OBJECTS;
      return OBJT(theElement)==BEOBJ && SIDE_ON_BND(theElement, i);
    }

    //! Returns true if at least one face of the element is a boundary face
    static bool isBoundaryElement(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::BEOBJ;
      using UG_NAMESPACE ::GM_OBJECTS;
      using UG::UINT;
      return OBJT(theElement)==BEOBJ;
    }

    //! \todo Please doc me!
    static int Edges_Of_Elem(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return EDGES_OF_ELEM(theElement);
    }

    //! \todo Please doc me!
    static int Corners_Of_Elem(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNERS_OF_ELEM(theElement);
    }

    //! \todo Please doc me!
    // Dummy implementation for vertices
    static int Corners_Of_Elem(const UG_NS< UG_DIM >::Node* theNode) {
      return 1;
    }

    //! Return number of corners of a given side
    static int Corners_Of_Side(const UG_NS< UG_DIM >::Element* theElement, int side) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNERS_OF_SIDE(theElement, side);
    }

    //! Return local number of a given corner of a given element side
    static int Corner_Of_Side(const UG_NS< UG_DIM >::Element* theElement, int side, int corner) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNER_OF_SIDE(theElement, side, corner);
    }

    //! Return local number of a given corner of a given element edge
    static int Corner_Of_Edge(const UG_NS< UG_DIM >::Element* theElement, int edge, int corner) {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      return CORNER_OF_EDGE(theElement, edge, corner);
    }

    //! Return number of sons of an element
    static int nSons(const UG_NAMESPACE ::element* element) {
      return UG_NAMESPACE ::ReadCW(element, UG_NAMESPACE ::NSONS_CE);
    }

    static int GetSons(const UG_NAMESPACE ::element* element, UG_NAMESPACE ::element* sonList[MAX_SONS]) {
      return UG_NAMESPACE ::GetSons(element, sonList);
    }

    /** \todo Remove the const casts */
    static int GetNodeContext(const UG_NAMESPACE ::element* element, const UG_NAMESPACE ::node** context) {
      return UG_NAMESPACE ::GetNodeContext(element, const_cast<UG_NAMESPACE ::node**>(context));
    }

    //! Encapsulates the GRID_ATTR macro
    static int Grid_Attr(const UG_NS< UG_DIM >::Grid* grid) {
      return GRID_ATTR(grid);
    }

    static int MarkForRefinement(UG_NAMESPACE ::element* element, int rule, int data) {
      return UG_NAMESPACE ::MarkForRefinement(element, (UG_NAMESPACE ::RefinementRule)rule, data);
    }

    //! Encapsulates the TAG macro
    static unsigned int Tag(const UG_NS< UG_DIM >::Element* theElement) {
      using UG::UINT;
      return TAG(theElement);
    }

    //! Doesn't ever get called, but needs to be there to calm the compiler
    static unsigned int Tag(const UG_NS< UG_DIM >::Node* theNode) {
      DUNE_THROW(GridError, "Called method Tag() for a vertex.  This should never happen!");
      return 0;
    }

    //! get corner in local coordinates, corner number in UG's numbering system
    template<class T>
    static void  getCornerLocal (const UG_NS< UG_DIM >::Element* theElement, int corner, FieldVector<T, UG_DIM>& local)
    {
      using UG_NAMESPACE ::element_descriptors;
      using UG::UINT;
      for (int i=0; i<UG_DIM; i++)
        local[i] = LOCAL_COORD_OF_TAG(TAG(theElement),corner)[i];
    }

    //! Next element in the UG element lists
    static UG_NS< UG_DIM >::Element* succ(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.succ;
    }

    //! Next element in the UG nodes lists
    static UG_NS< UG_DIM >::Node* succ(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->succ;
    }

    //! Calm the compiler
    static void* succ(const void* theWhatever) {
      DUNE_THROW(NotImplemented, "No successor available for this kind of object");
      return 0;
    }

    //! Return true if the element is a leaf element
    static bool isLeaf(const UG_NS< UG_DIM >::Element* theElement) {
      return UG_NAMESPACE ::EstimateHere(theElement);
    }

    //! Return true if the node is a leaf node
    static bool isLeaf(const UG_NS< UG_DIM >::Node* theNode) {
#ifndef ModelP
#warning Method isLeaf() for nodes will not work properly in case of vertical load balancing
#endif
      return !theNode->son;
    }

    // /////////////////////////////////////////////
    //   Level indices
    // /////////////////////////////////////////////

    //! Gets the level index of a UG element
    static int& levelIndex(UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.levelIndex;
    }

    //! Gets the level index of a UG element
    static const int& levelIndex(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.levelIndex;
    }

    //! Gets the level index of a UG sidevector
    static int& levelIndex(Vector* theVector) {
#if UG_DIM == 2
      DUNE_THROW(GridError, "levelIndex in side vector only in 3D!");
#endif
      return reinterpret_cast<int&>(theVector->index);
    }

    //! Gets the level index of a UG sidevector
    static const int& levelIndex(const Vector* theVector) {
#if UG_DIM == 2
      DUNE_THROW(GridError, "levelIndex in side vector only in 3D!");
#endif
      return reinterpret_cast<const int&>(theVector->index);
    }

    //! Gets the level index of a UG edge
    static int& levelIndex(UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->levelIndex;
    }

    //! Gets the level index of a UG edge
    static const int& levelIndex(const UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->levelIndex;
    }

    //! Gets the level index of a UG node
    static int& levelIndex(UG_NS< UG_DIM >::Node* theNode) {
      return theNode->levelIndex;
    }

    //! Gets the level index of a UG node
    static const int& levelIndex(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->levelIndex;
    }

    // /////////////////////////////////////////////
    //   Leaf indices
    // /////////////////////////////////////////////

    //! Gets the leaf index of a UG element
    static int& leafIndex(UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.leafIndex;
    }

    //! Gets the leaf index of a UG element
    static const int& leafIndex(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.leafIndex;
    }

    //! Gets the leaf index of a UG sidevector
    static int& leafIndex(Vector* theVector) {
      return reinterpret_cast<int &>(theVector->skip);
    }

    //! Gets the leaf index of a UG sidevector
    static const int& leafIndex(const Vector* theVector) {
      return reinterpret_cast<const int &>(theVector->skip);
    }

    //! Gets the leaf index of a UG edge
    static int& leafIndex(UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->leafIndex;
    }

    //! Gets the leaf index of a UG edge
    static const int& leafIndex(const UG_NS< UG_DIM >::Edge* theEdge) {
      return theEdge->leafIndex;
    }

    //! Gets the leaf index of a UG node
    static int& leafIndex(UG_NS< UG_DIM >::Node* theNode) {
      return theNode->myvertex->iv.leafIndex;
    }

    //! Gets the leaf index of a UG node
    static const int& leafIndex(const UG_NS< UG_DIM >::Node* theNode) {
      return theNode->myvertex->iv.leafIndex;
    }

    // /////////////////////////////////////////////
    //   IDs
    // /////////////////////////////////////////////

    //! Gets the index of a UG element
    static unsigned int id(const UG_NS< UG_DIM >::Element* theElement) {
      return theElement->ge.id;
    }

    //! Gets the index of a UG node
    static unsigned int id(const UG_NS< UG_DIM >::Node* theNode) {
#if UG_DIM == 2
      return theNode->myvertex->iv.id | 0x80000000;
#else
      return theNode->myvertex->iv.id | 0xC0000000;
#endif
    }

    //! \todo Please doc me!
    static void Local_To_Global(int n, double** y,
                                const FieldVector<double, UG_DIM>& local,
                                FieldVector<double, UG_DIM>& global) {
      using UG::DOUBLE;
      LOCAL_TO_GLOBAL(n,y,local,global);
    }

    /**
     * \param n Number of corners of the element
     * \param x Coordinates of the corners of the element
     * \param local Local evaluation point
     *
     * \return The return type is int because the macro INVERSE_TRANSFORMATION
     *  returns 1 on failure.
     */
    static int Transformation(int n, double** x,
                              const FieldVector<double, UG_DIM>& local, FieldMatrix<double,UG_DIM,UG_DIM>& mat) {
      using UG_NAMESPACE ::DOUBLE_VECTOR;
      using UG::DOUBLE;
      double det;
#ifndef SMALL_D
      const double SMALL_D = DBL_EPSILON*10;
#endif
      INVERSE_TRANSFORMATION(n, x, local, mat, det);
      return 0;
    }

    /** \brief Compute the integration element on an element face
        \param nc Number of corners of the element face
        \param co_global Coordinates of the corners of the face
        \param ip_local Local coordinates where integration element is to be evaluated
     */
    static double SurfaceElement(int nc,
                                 const double co_global[MAX_CORNERS_OF_ELEM][ UG_DIM ],
                                 const double ip_local[ UG_DIM ]) {
      double result;
      if (UG_NAMESPACE ::SurfaceElement(UG_DIM, nc, co_global, ip_local, &result))
        DUNE_THROW(GridError, "UG_D" << UG_DIM << "::SurfaceElement returned error code!");
      return result;
    }

    //! Returns the i-th corner of a UG element
    static UG_NS< UG_DIM >::Node* Corner(const UG_NS< UG_DIM >::Element* theElement, int i) {
      using UG_NAMESPACE ::NODE;
      using UG_NAMESPACE ::n_offset;
      using UG::UINT;
      return CORNER(theElement, i);
    }

    //! get edge from node i to node j (in UG's numbering !
    static UG_NS< UG_DIM >::Edge* GetEdge (UG_NS< UG_DIM >::Node* nodei, UG_NS< UG_DIM >::Node* nodej) {
      return UG_NAMESPACE ::GetEdge(nodei,nodej);
    }

    //! access side vector from element (this is just a dummy to compile code also in 2d)
    static Vector* SideVector (const UG_NS< UG_DIM >::Element* theElement, int i)
    {
#if UG_DIM == 2
      DUNE_THROW(GridError, "side vector only in 3D!");
#else
      using UG::D3::VECTOR;
      using UG::D3::svector_offset;
      using UG::UINT;
      return SVECTOR(theElement,i);
#endif
    }

    /** \brief Return a pointer to the father of the given element */
    static UG_NS< UG_DIM >::Element* EFather(const UG_NS< UG_DIM >::Element* theElement) {
      using UG_NAMESPACE ::ELEMENT;
      using UG_NAMESPACE ::father_offset;
      using UG::UINT;
      return EFATHER(theElement);
    }

    //! get father element of vertex
    static UG_NS< UG_DIM >::Element* NFather(UG_NS< UG_DIM >::Node* theNode) {
      return theNode->myvertex->iv.father;
    }

    //! get father node of vertex
    static UG_NS< UG_DIM >::Node* NodeNodeFather(UG_NS< UG_DIM >::Node* theNode) {
      using UG_NAMESPACE ::NDOBJ;
      using UG_NAMESPACE ::GM_OBJECTS;
      using UG::UINT;
      if (theNode->father==0)
        return 0;         // no father at all
      if (OBJT(theNode->father)==NDOBJ)
        return (UG_NAMESPACE ::node*) theNode->father;
      else
        return 0;         // may be edge or element
    }

    static unsigned int ReadCW(void* obj, int ce) {
      return UG_NAMESPACE ::ReadCW(obj, ce);
    }

    static void WriteCW(void* obj, int ce, int n) {
      UG_NAMESPACE ::WriteCW(obj, ce, n);
    }

    //! \todo Please doc me!
    static int InitUg(int* argcp, char*** argvp) {
      return UG_NAMESPACE ::InitUg(argcp, argvp);
    }

    static void ExitUg() {
      UG_NAMESPACE ::ExitUg();
    }

    static void DisposeMultiGrid(UG_NAMESPACE ::multigrid* mg) {
      UG_NAMESPACE ::DisposeMultiGrid(mg);
    }

#ifndef UG_LGMDOMAIN
    //! \todo Please doc me!
    static void* CreateBoundaryValueProblem(const char* BVPname,
                                            int numOfCoeffFunc,
                                            UG_NAMESPACE ::CoeffProcPtr coeffs[],
                                            int numOfUserFct,
                                            UG_NAMESPACE ::UserProcPtr userfct[]) {
      return UG_NAMESPACE ::CreateBoundaryValueProblem(BVPname, 0, numOfCoeffFunc, coeffs,
                                                       numOfUserFct, userfct);
    }
#endif

    //! Set the current boundary value problem
    static void Set_Current_BVP(void** thisBVP) {
      UG_NAMESPACE ::Set_Current_BVP(thisBVP);
    }

    //! Get UG boundary value problem from its name
    static void** BVP_GetByName(const char* name) {
      return UG_NAMESPACE ::BVP_GetByName(name);
    }

    //! Dispose of a boundary value problem
    static int BVP_Dispose(void** BVP) {
      return UG_NAMESPACE ::BVP_Dispose(BVP);
    }

    //! Get UG multigrid object from its name
    static UG_NS< UG_DIM >::MultiGrid* GetMultigrid(const char* name) {
      return UG_NAMESPACE ::GetMultigrid(name);
    }

    static int LBCommand(int argc, const char** argv) {
      /** \todo Can we remove the cast? */
      return UG_NAMESPACE ::LBCommand(argc, (char**)argv);
    }

    static int ConfigureCommand(int argc, const char** argv) {
      /** \todo Kann man ConfigureCommand so �ndern da� man auch ohne den cast auskommt? */
      return UG_NAMESPACE ::ConfigureCommand(argc, (char**)argv);
    }

    static int NewCommand(int argc, char** argv) {
      return UG_NAMESPACE ::NewCommand(argc, argv);
    }

    static int CreateFormatCmd(int argc, char** argv) {
      return UG_NAMESPACE ::CreateFormatCmd(argc, argv);
    }

#ifndef UG_LGMDOMAIN
    static void* CreateDomain(const char* name, const double* midPoint, double radius,
                              int segments, int corners, int convex) {
      return UG_NAMESPACE ::CreateDomain(name, midPoint, radius, segments, corners, convex);
    }
#endif

    static void* InsertInnerNode(UG_NAMESPACE ::grid* grid, const double* pos) {
      return UG_NAMESPACE ::InsertInnerNode(grid, pos);
    }

#ifndef UG_LGMDOMAIN
    static void* CreateBoundarySegment(const char *name, int left, int right,
                                       int index, int res,
                                       int *point,
                                       const double *alpha, const double *beta,
                                       UG_NAMESPACE ::BndSegFuncPtr boundarySegmentFunction,
                                       void *userData) {
      return UG_NAMESPACE ::CreateBoundarySegment(name,            // internal name of the boundary segment
                                                  left,             //  id of left subdomain
                                                  right,             //  id of right subdomain
                                                  index,         // Index of the segment
                                                  UG_NAMESPACE ::NON_PERIODIC, // I don't know what this means
                                                  res,             // Resolution, only for the UG graphics
                                                  point,
                                                  alpha,
                                                  beta,
                                                  boundarySegmentFunction,
                                                  userData);
    }
#endif

  };

  template <>
  class UG_NS< UG_DIM >::Entity<0>
  {
  public:
    typedef UG_NAMESPACE ::element T;
  };

  template <>
  class UG_NS< UG_DIM >::Entity< UG_DIM - 1 >
  {
  public:
    typedef UG_NAMESPACE ::edge T;
  };

  template <>
  class UG_NS< UG_DIM >::Entity< UG_DIM > {
  public:
    typedef UG_NAMESPACE ::node T;
  };

#undef UG_NAMESPACE

#if defined ModelP && UG_DIM==3
#undef PRIO2LISTPART
#endif

} // namespace Dune
