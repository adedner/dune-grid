// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <iostream>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/test/checkvtkfile.hh>

// Test FieldInfo::writeComponents static method
int testWriteComponents()
{
  using FI = Dune::VTK::FieldInfo;
  using Type = FI::Type;
  int result = 0;

  // unspecified: returns size as-is
  if (FI::writeComponents(Type::unspecified, 1) != 1) result = 1;
  if (FI::writeComponents(Type::unspecified, 2) != 2) result = 1;
  if (FI::writeComponents(Type::unspecified, 4) != 4) result = 1;
  if (FI::writeComponents(Type::unspecified, 5) != 5) result = 1;

  // scalar: always returns 1
  if (FI::writeComponents(Type::scalar, 1) != 1) result = 1;

  // vector: always returns 3, throws if size > 3
  if (FI::writeComponents(Type::vector, 1) != 3) result = 1;
  if (FI::writeComponents(Type::vector, 2) != 3) result = 1;
  if (FI::writeComponents(Type::vector, 3) != 3) result = 1;

  // vector with size > 3 should throw
  try {
    FI::writeComponents(Type::vector, 4);
    result = 1; // should not reach here
  } catch (const Dune::IOError&) {
    // expected
  }

  // tensor: should throw NotImplemented
  try {
    FI::writeComponents(Type::tensor, 9);
    result = 1; // should not reach here
  } catch (const Dune::NotImplemented&) {
    // expected
  }

  return result;
}

template< class GridView >
class VTK2DVectorFunction
  : public Dune :: VTKWriter< GridView > :: VTKFunction
{
  typedef typename GridView :: Grid :: ctype DT;
  typedef typename GridView :: template Codim< 0 > :: Entity Entity;
public:
  VTK2DVectorFunction() {}

  virtual int ncomps () const { return 2; }

  virtual double evaluate (int comp, [[maybe_unused]] const Entity& e,
                           [[maybe_unused]] const Dune::FieldVector<DT,GridView::dimension>& xi) const
  {
    return comp * 1.0;
  }

  virtual std::string name () const
  {
    return "2D-vector";
  }
};

template< class GridView >
int testFieldInfoTypes(Dune::VTKChecker& vtkChecker, const GridView& gridView)
{
  constexpr static int dim = GridView :: dimension;

  Dune :: VTKWriter< GridView > vtk( gridView, Dune::VTK::conforming );

  // Test unspecified with 4 components
  std::size_t numVertices = gridView.indexSet().size(dim);
  std::vector<double> unspecifiedData4(numVertices * 4, 0.0);
  for (std::size_t i = 0; i < unspecifiedData4.size(); ++i)
    unspecifiedData4[i] = static_cast<double>(i % 4);
  vtk.addVertexData(unspecifiedData4, "unspecified-4comp", 4);

  // Test scalar
  auto scalarFunc = [](const Dune::FieldVector<double,dim>& x) { return x.two_norm(); };
  vtk.addCellData(scalarFunc,
                  Dune::VTK::FieldInfo("scalar-lambda",
                                       Dune::VTK::FieldInfo::Type::scalar, 1));

  // Test vector with 2D data (should be padded to 3 in output)
  vtk.addVertexData(std::make_shared<VTK2DVectorFunction<GridView>>());

  // Test vector with dim components
  auto vecFunc = [](const Dune::FieldVector<double,dim>& x) { return x; };
  vtk.addVertexData(vecFunc,
                    Dune::VTK::FieldInfo("vector-lambda",
                                         Dune::VTK::FieldInfo::Type::vector, dim));

  std::string name;
  std::ostringstream prefix;
  prefix << "fieldinfo-test-" << dim << "D";
  int rank = gridView.comm().rank();

  name = vtk.write(prefix.str() + "-ascii");
  if(rank == 0) vtkChecker.push(name);

  name = vtk.write(prefix.str() + "-base64", Dune::VTK::base64);
  if(rank == 0) vtkChecker.push(name);

  return 0;
}

int main(int argc, char **argv)
{
  try {
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if(mpiHelper.rank() == 0)
      std::cout << "fieldinfo-test: MPI_Comm_size == " << mpiHelper.size() << std::endl;

    int result = 0;

    // Test writeComponents static method
    result = testWriteComponents();
    if (result != 0) {
      std::cerr << "testWriteComponents failed!" << std::endl;
      return result;
    }
    std::cout << "testWriteComponents passed." << std::endl;

    Dune::VTKChecker vtkChecker;

    // Test FieldInfo types with actual VTK file output
    {
      Dune::YaspGrid<2> g2({1.0, 2.0}, {4, 4});
      g2.globalRefine(1);
      testFieldInfoTypes(vtkChecker, g2.leafGridView());
    }
    {
      Dune::YaspGrid<3> g3({1.0, 2.0, 3.0}, {4, 4, 4});
      g3.globalRefine(1);
      testFieldInfoTypes(vtkChecker, g3.leafGridView());
    }

    result = vtkChecker.check();

    return result;

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }
}
