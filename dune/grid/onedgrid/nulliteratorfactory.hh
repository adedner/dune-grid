// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ONEDGRID_NULL_ITERATORS_HH
#define DUNE_ONEDGRID_NULL_ITERATORS_HH

#include "onedgridlist.hh"

namespace Dune {

  template <int mydim, int cdim, class ct> class OneDEntityImp;

  template <int dim, int dimw, class ct>
  class OneDGridNullIteratorFactory {};

  template <int cdim, class ct>
  class OneDGridNullIteratorFactory<0,cdim,ct> {

  public:

    static typename OneDGridList<OneDEntityImp<0,cdim,ct> >::iterator null() {
      return emptyList_.end();
    }

  private:
    static OneDGridList<OneDEntityImp<0,cdim,ct> > emptyList_;
  };

  template <int cdim, class ct>
  class OneDGridNullIteratorFactory<1,cdim,ct> {

  public:

    static typename OneDGridList<OneDEntityImp<1,cdim,ct> >::iterator null() {
      return emptyList_.end();
    }

  private:
    static OneDGridList<OneDEntityImp<1,cdim,ct> > emptyList_;
  };

  extern template class OneDGridNullIteratorFactory<0,1,double>;
  extern template class OneDGridNullIteratorFactory<1,1,double>;
  extern template class OneDGridNullIteratorFactory<0,2,double>;
  extern template class OneDGridNullIteratorFactory<1,2,double>;
  extern template class OneDGridNullIteratorFactory<0,3,double>;
  extern template class OneDGridNullIteratorFactory<1,3,double>;

} // end namespace Dune

#endif
