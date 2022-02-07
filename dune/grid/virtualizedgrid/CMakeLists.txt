set(HEADERS
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
