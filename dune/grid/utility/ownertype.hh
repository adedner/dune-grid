// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_UTILITY_OWNERTYPE_HH
#define DUNE_GRID_UTILITY_OWNERTYPE_HH

namespace Dune {

/** \brief Attributes used in the generic ownership model */
enum class OwnerType {
  Owner=0,          //!< all interior and owned border entities
  Overlap=1,        //!< all entities lying in the overlap zone and non-owned border entities
  Ghost=2,          //!< all ghost and front entities
  Unknown=10        //!< No unique owner type can be specified
};

} // end namespace Dune

#endif // DUNE_GRID_UTILITY_OWNERTYPE_HH
