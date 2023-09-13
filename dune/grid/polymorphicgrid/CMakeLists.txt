# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  datahandle.hh
  entity.hh
  entityiterator.hh
  entityseed.hh
  geometry.hh
  grid.hh
  idtype.hh
  indexsets.hh
  intersectioniterator.hh
  intersections.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/polymorphicgrid/)

add_subdirectory(common)