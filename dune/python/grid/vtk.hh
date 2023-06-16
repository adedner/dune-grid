// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYTHON_GRID_VTK_HH
#define DUNE_PYTHON_GRID_VTK_HH

#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/common/getdimension.hh>
#include <dune/python/grid/enums.hh>
#include <dune/python/grid/object.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Python
  {

    // External Forward Declarations
    // -----------------------------

    template< class GridFunction >
    struct GridFunctionTraits;



    // addToVTKWriter
    // --------------

    template< class GridFunction, class Writer = VTKWriter< std::decay_t< decltype( gridView( std::declval< const GridFunction & >() ) ) > > >
    void addToVTKWriter ( const GridFunction &gf, std::string name, Writer &vtkWriter, VTKDataType dataType )
    {
      typedef typename GridFunctionTraits< GridFunction >::Range Range;

      VTK::FieldInfo::Type type = VTK::FieldInfo::Type::none;
      switch( dataType ) {
        case VTKDataType::CellScalar:
        case VTKDataType::PointScalar:
          type = VTK::FieldInfo::Type::scalar; break;
        case VTKDataType::CellVector:
        case VTKDataType::PointVector:
          type = VTK::FieldInfo::Type::vector; break;
        case VTKDataType::CellTensor:
        case VTKDataType::PointTensor:
          type = VTK::FieldInfo::Type::tensor; break;
        default:
          type = VTK::FieldInfo::Type::none;
      }

      VTK::FieldInfo info( std::move( name ), type, GetDimension<Range>::value );

      switch( dataType ) {
        case VTKDataType::CellData:
        case VTKDataType::CellScalar:
        case VTKDataType::CellVector:
        case VTKDataType::CellTensor:
          vtkWriter.addCellData( gf, info );
          break;
        default:
          vtkWriter.addVertexData( gf, info );
      }
    }



    // registerVTKWriter
    // -----------------

    template<class GridView>
    void registerVTKWriter(pybind11::handle scope)
    {
      {
        typedef VTKWriter< GridView > Writer;

        auto cls = insertClass<Writer>(scope, "VTKWriter",
            GenerateTypeName("VTKWriter", MetaType<GridView>()) ).first;

        cls.def( "write",
            [] ( Writer &writer, const std::string &name, Dune::VTK::OutputType outputType ) {
              writer.write( name, outputType );
            },
            pybind11::arg("name"),
            pybind11::arg("outputType")=VTK::appendedraw );

        cls.def( "write",
            [] ( Writer &writer, const std::string &name, int number, Dune::VTK::OutputType outputType ) {
              std::stringstream s; s << name << std::setw(5) << std::setfill('0') << number;
              writer.write( s.str(), outputType );
            },
            pybind11::arg("name"),
            pybind11::arg("number"),
            pybind11::arg("outputType")=VTK::appendedraw );
      }

      {
        typedef SubsamplingVTKWriter< GridView > Writer;

        auto cls = insertClass< Writer, VTKWriter<GridView> >(scope, "SubsamplingVTKWriter",
            GenerateTypeName("SubsamplingVTKWriter", MetaType<GridView>()) ).first;

        cls.def( "write",
            [] ( Writer &writer, const std::string &name, Dune::VTK::OutputType outputType ) {
              writer.write( name, outputType );
            },
            pybind11::arg("name"),
            pybind11::arg("outputType")=VTK::appendedraw );

        cls.def( "write",
            [] ( Writer &writer, const std::string &name, int number ) {
              std::stringstream s; s << name << std::setw(5) << std::setfill('0') << number;
              writer.write( s.str() );
            },
            pybind11::arg("name"),
            pybind11::arg("number") );
        cls.def( "write",
            [] ( Writer &writer, const std::string &name, int number, Dune::VTK::OutputType outputType ) {
              std::stringstream s; s << name << std::setw(5) << std::setfill('0') << number;
              writer.write( s.str(), outputType );
            },
            pybind11::arg("name"),
            pybind11::arg("number"),
            pybind11::arg("outputType")=VTK::appendedraw );
      }
    }

  } // namespace Python

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_VTK_HH
