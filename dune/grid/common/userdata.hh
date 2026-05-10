// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_COMMON_USERDATA_HH
#define DUNE_GRID_COMMON_USERDATA_HH

#include <typeinfo>
#include <memory>

namespace Dune
{
  template <class Implementation>
  class GridUserDataObject;

  /** \brief interface class for user data storage
   */
  class GridUserDataIF
  {
  protected:
    GridUserDataIF () {}
  public:
    virtual ~GridUserDataIF() {}
    virtual std::string name() const { return std::string("interface"); }

    template <class Implementation>
    bool compare() const
    {
      return typeid(Implementation).name() == this->name();
    }

    template <class Implementation>
    std::shared_ptr< Implementation >& get()
    {
      typedef GridUserDataObject< Implementation > GridUserDataObjectType;
      //assert( dynamic_cast< GridUserDataObjectType& > (*this) );
      return dynamic_cast< GridUserDataObjectType& > (*this).object();
    }
  };

  template <class Implementation>
  class GridUserDataObject : public GridUserDataIF
  {
    typedef std::shared_ptr< Implementation > ObjectPointerType;
    ObjectPointerType obj_;
  public:
    GridUserDataObject( const ObjectPointerType& obj )
      : obj_( obj )
    {}

    GridUserDataObject( ObjectPointerType&& obj )
      : obj_( std::forward(obj) )
    {}

    std::string name() const override { return typeid(Implementation).name(); }

    std::shared_ptr< Implementation >& object() { return obj_; }
  };

} // end namespace Dune

#endif
