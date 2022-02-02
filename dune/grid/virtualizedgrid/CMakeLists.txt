set(HEADERS
  virtualizedgrvirtualized.hh
  virtualizedgrvirtualizedseed.hh
  virtualizedgridhierarchiciterator.hh
  virtualizedgridgeometry.hh
  virtualizedgridintersectioniterator.hh
  virtualizedgridintersections.hh
  virtualizedgridindexsets.hh
  virtualizedgridleafiterator.hh
  virtualizedgridleveliterator.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/grid/virtualizedgrid)
