// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "nulliteratorfactory.hh"

namespace Dune {

template<int dimw, class ct>
OneDGridList<OneDEntityImp<0,dimw,ct> > OneDGridNullIteratorFactory<0,dimw,ct>::emptyList_;
template<int dimw, class ct>
OneDGridList<OneDEntityImp<1,dimw,ct> > OneDGridNullIteratorFactory<1,dimw,ct>::emptyList_;

template class OneDGridNullIteratorFactory<0,1,double>;
template class OneDGridNullIteratorFactory<1,1,double>;
template class OneDGridNullIteratorFactory<0,2,double>;
template class OneDGridNullIteratorFactory<1,2,double>;
template class OneDGridNullIteratorFactory<0,3,double>;
template class OneDGridNullIteratorFactory<1,3,double>;

} // end namespace Dune
