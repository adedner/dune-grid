# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

set(HEADERS
  cast.hh
  datahandle.hh
  entity.hh
  entityseed.hh
  geometry.hh
  hierarchiciterator.hh
  idtype.hh
  indexsets.hh
  intersectioniterator.hh
  intersections.hh
  leafiterator.hh
  leveliterator.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/virtualizedgrid/)
