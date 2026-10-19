// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_TEST_CHECKCOMMCORRECTNESS_HH
#define DUNE_GRID_TEST_CHECKCOMMCORRECTNESS_HH

#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>

#include <dune/common/hash.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/test/coordinatetolerance.hh>
#include <dune/grid/test/periodicmap.hh>

namespace Dune {

  std::size_t hash_value(const PartitionType& pt)
  {
    return static_cast<std::size_t>(pt);
  }

}

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(),DUNE_HASH_TYPE(Dune::PartitionType))

namespace Dune {

  namespace GridCheck {

      struct SymmetryVerifyingDataHandle
        : public Dune::CommDataHandleIF<SymmetryVerifyingDataHandle,std::size_t>
      {

        bool contains(int /* dim */, int codim) const
        {
          return codim == _codim;
        }

        bool fixedSize(int /* dim */, int /* codim */) const
        {
          return _fixed_size;
        }

        template<typename E>
        std::size_t size(const E&) const
        {
          return _writes_per_cell;
        }

        template<typename Buf, typename E>
        void gather(Buf& buf, const E&) const
        {
          for (std::size_t i = 0; i < _writes_per_cell; ++i)
            buf.write(std::size_t(42));
          _writes += _writes_per_cell;
        }

        template<typename Buf, typename E>
        void scatter(Buf& buf, const E&, std::size_t n) const
        {
          assert(_writes_per_cell == n);
          for (std::size_t i = 0; i < _writes_per_cell; ++i)
            {
              std::size_t tmp = 0;
              buf.read(tmp);
              assert(tmp == 42);
            }
          _reads += _writes_per_cell;
        }

        SymmetryVerifyingDataHandle(int codim, bool fixed_size, std::size_t writes_per_cell)
          : _codim(codim)
          , _fixed_size(fixed_size)
          , _writes_per_cell(writes_per_cell)
          , _reads(0)
          , _writes(0)
        {}

        const int _codim;
        const bool _fixed_size;
        const std::size_t _writes_per_cell;
        mutable std::size_t _reads;
        mutable std::size_t _writes;

      };

      struct CodimLayout
      {

        bool operator()(Dune::GeometryType gt, int dim) const
        {
          return (static_cast< int >( gt.dim() ) == dim - _codim);
        }

        int _codim;

      };

      /** \brief Verify that communicating entities have the same global id
       *
       * The sender writes the hash of the global id of the entity, the receiver
       * compares it with the hash of the global id of its own entity. In
       * parallel this includes periodic copies whose original is stored on
       * another process, which the local checks cannot see.
       */
      template<typename GV>
      struct IdVerifyingDataHandle
        : public Dune::CommDataHandleIF<IdVerifyingDataHandle<GV>, std::size_t>
      {
        using IdType = typename GV::Grid::GlobalIdSet::IdType;

        bool contains(int /* dim */, int codim) const
        {
          return codim == _codim;
        }

        bool fixedSize(int /* dim */, int /* codim */) const
        {
          return true;
        }

        template<typename E>
        std::size_t size(const E&) const
        {
          return 1;
        }

        template<typename Buf, typename E>
        void gather(Buf& buf, const E& e) const
        {
          buf.write(std::hash<IdType>{}(_gv.grid().globalIdSet().id(e)));
        }

        template<typename Buf, typename E>
        void scatter(Buf& buf, const E& e, std::size_t /* n */)
        {
          std::size_t hash = 0;
          buf.read(hash);
          if (hash != std::hash<IdType>{}(_gv.grid().globalIdSet().id(e)))
          {
            std::cerr << "[" << _gv.comm().rank() << "] Error: received a different global id for codim "
                      << _codim << " entity at " << e.geometry().center() << std::endl;
            ++_errors;
          }
        }

        IdVerifyingDataHandle(GV gv, int codim)
          : _gv(gv), _codim(codim)
        {}

        GV _gv;
        const int _codim;
        std::size_t _errors = 0;
      };

      template<typename GV>
      struct CommunicationTestDataHandle
        : public Dune::CommDataHandleIF<CommunicationTestDataHandle<GV>,
                                  typename GV::template Codim<0>::Geometry::GlobalCoordinate>
      {
        using ctype = typename GV::ctype;

        bool contains([[maybe_unused]] int dim, int codim) const
        {
          return codim == _codim;
        }

        bool fixedSize([[maybe_unused]] int dim, [[maybe_unused]] int codim) const
        {
          return true;
        }

        template<typename E>
        std::size_t size([[maybe_unused]] const E& e) const
        {
          return 1;
        }

        // debug output is collected in _log and only printed if a check fails
        void require(bool condition, const std::string& what) const
        {
          if (!condition)
            std::cerr << "[" << _gv.comm().rank() << "] Error in communication test: " << what << "\n"
                      << "[" << _gv.comm().rank() << "] debug output of this communication:\n"
                      << _log.str() << std::flush;
          assert(condition);
        }

        // send the center of the entity, mapped to the periodic unit cell
        template<typename Buf, typename E>
        void gather(Buf& buf, const E& e) const
        {
          auto center = e.geometry().center();
          _log << "Gathering from entity " << _mapper.index(e) << " at " << center << "\n";
          require(_allowed_writes.count(e.partitionType()) > 0,
                  "write from a " + PartitionName(e.partitionType()) + " entity is not allowed");
          buf.write(_periodicMap(center));
          ++_writes[_mapper.index(e)];
          require(GridCheck::coordinatesEqual(center, _coords[_mapper.index(e)], _tolerance),
                  "center of the gathered entity differs from the stored one");
        }

        // the received center must match the own center, up to the periodic map
        template<typename Buf, typename E>
        void scatter(Buf& buf, const E& e, [[maybe_unused]] std::size_t n) const
        {
          auto center = e.geometry().center();
          _log << "Scattering to entity " << _mapper.index(e) << " at " << center << "\n";
          require(_allowed_reads.count(e.partitionType()) > 0,
                  "read into a " + PartitionName(e.partitionType()) + " entity is not allowed");
          typename E::Geometry::GlobalCoordinate data;
          buf.read(data);
          ++_reads[_mapper.index(e)];
          require(GridCheck::coordinatesEqual(data, _periodicMap(center), _tolerance),
                  "received center differs from the own center (up to the periodic map)");
          require(GridCheck::coordinatesEqual(center, _coords[_mapper.index(e)], _tolerance),
                  "center of the scattered entity differs from the stored one");
        }

        // check the numbers of writes and reads of each entity; check returns whether they are correct
        template<int cd, typename Check>
        void verify(Dune::Codim<cd> codim, Check check)
        {
          assert(codim == _codim);
          for (const auto& e : entities(_gv,codim))
            {
              std::ostringstream entity;
              entity << "Entity of codim " << cd
                     << " at " << e.geometry().center()
                     << " with partition " << PartitionName(e.partitionType())
                     << " and " << _writes[_mapper.index(e)] << " / " << _reads[_mapper.index(e)] << " writes / reads";
              _log << entity.str() << "\n";
              require(check(
                        e.partitionType(),
                        _allowed_writes.count(e.partitionType()) > 0,
                        _allowed_reads.count(e.partitionType()) > 0,
                        _writes[_mapper.index(e)],
                        _reads[_mapper.index(e)]
                        ),
                      "unexpected number of writes / reads: " + entity.str());
            }
        }

        CommunicationTestDataHandle(GV gv, int codim, const std::unordered_set<Dune::PartitionType>& allowed_writes, const std::unordered_set<Dune::PartitionType>& allowed_reads, const std::vector<typename GV::template Codim<0>::Geometry::GlobalCoordinate>& coords,
                                    const PeriodicMap<GV>& periodicMap = noPeriodicMap<GV>())
          : _gv(gv)
          , _codim(codim)
          , _mapper(gv, CodimLayout{_codim})
          , _allowed_writes(allowed_writes)
          , _allowed_reads(allowed_reads)
          , _reads(_mapper.size(),0)
          , _writes(_mapper.size(),0)
          , _coords(coords)
          , _tolerance(GridCheck::coordinateTolerance(gv))
          , _periodicMap(periodicMap)
        {}

        GV _gv;
        const int _codim;
        Dune::MultipleCodimMultipleGeomTypeMapper<GV> _mapper;
        const std::unordered_set<Dune::PartitionType> _allowed_writes;
        const std::unordered_set<Dune::PartitionType> _allowed_reads;
        mutable std::vector<std::size_t> _reads;
        mutable std::vector<std::size_t> _writes;
        const std::vector<typename GV::template Codim<0>::Geometry::GlobalCoordinate>& _coords;
        const ctype _tolerance;
        const PeriodicMap<GV> _periodicMap;
        mutable std::ostringstream _log;

      };

      template<typename GV, int cd>
      void check_communication_correctness_do(GV gv, Codim<cd> codim, const PeriodicMap<GV>& periodicMap)
      {
        if (gv.grid().comm().rank() == 0)
          {
            std::cout << "Checking codim " << cd << std::endl;
          }

        std::unordered_map<Dune::PartitionType,std::size_t> count;

        Dune::MultipleCodimMultipleGeomTypeMapper<GV> mapper(gv, CodimLayout{codim});

        std::vector<
          typename GV::template Codim<0>::Geometry::GlobalCoordinate
          > coords(mapper.size());

        // whether this process stores several copies of an entity, i.e. is
        // its own neighbor across a periodic boundary (a copy whose original
        // is stored on another process does not count)
        const bool hasLocalCoincidences = PeriodicIdentification<GV>(gv, cd, periodicMap).hasCoincidences();

        // start by counting entities by partition type and storing entity positions
        for (const auto& e : entities(gv,codim))
          {
            ++count[e.partitionType()];
            coords[mapper.index(e)] = e.geometry().center();
          }

        // communicating entities must have the same global id; in parallel this
        // checks the periodic copies whose original is stored on another process
        {
          IdVerifyingDataHandle<GV> dh(gv, cd);
          gv.communicate(dh, InteriorBorder_All_Interface, ForwardCommunication);
          assert(dh._errors == 0);
          gv.communicate(dh, InteriorBorder_All_Interface, BackwardCommunication);
          assert(dh._errors == 0);
        }

        {
          SymmetryVerifyingDataHandle dh_forward(cd,false,3);
          gv.communicate(dh_forward,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);

          SymmetryVerifyingDataHandle dh_backward(cd,true,3);
          gv.communicate(dh_backward,InteriorBorder_InteriorBorder_Interface,BackwardCommunication);

          if (count[BorderEntity] > 0)
            {
              assert(dh_forward._writes == dh_forward._reads);
              assert(dh_backward._writes == dh_backward._reads);

              assert(dh_forward._writes == dh_backward._writes);
              assert(dh_backward._writes == dh_backward._writes);
            }

          // with two processes each border entity is shared with exactly one other
          // process, unless a process is its own neighbor across a periodic boundary
          if (gv.grid().comm().size() == 2 && !hasLocalCoincidences)
            {
              if (gv.comm().rank() == 0)
                std::cout << "MPI size == 2, checking writes (" << (dh_forward._writes / 3)
                          << ") against count of border entities (" <<  count[BorderEntity] << ")"
                          << std::endl;
              assert((dh_forward._writes / 3) == count[BorderEntity]);
            }

        }

        {
          SymmetryVerifyingDataHandle dh_forward(cd,true,7);
          gv.communicate(dh_forward,All_All_Interface,ForwardCommunication);

          SymmetryVerifyingDataHandle dh_backward(cd,false,7);
          gv.communicate(dh_backward,All_All_Interface,BackwardCommunication);

          if (count[BorderEntity] > 0)
            {
              assert(dh_forward._writes == dh_forward._reads);
              assert(dh_backward._writes == dh_backward._reads);

              assert(dh_forward._writes == dh_backward._writes);
              assert(dh_backward._writes == dh_backward._writes);
            }

        }

        using PTSet = std::unordered_set<PartitionType>;

        {

          PTSet writers({InteriorEntity,BorderEntity});
          PTSet readers({InteriorEntity,BorderEntity,OverlapEntity,FrontEntity,GhostEntity});

          CommunicationTestDataHandle<GV> dh(gv,codim,writers,readers,coords,periodicMap);
          gv.communicate(dh,InteriorBorder_All_Interface,ForwardCommunication);

        }

        {

          PTSet writers({InteriorEntity,BorderEntity,OverlapEntity,FrontEntity,GhostEntity});
          PTSet readers({InteriorEntity,BorderEntity});

          CommunicationTestDataHandle<GV> dh(gv,codim,writers,readers,coords,periodicMap);
          gv.communicate(dh,InteriorBorder_All_Interface,BackwardCommunication);

        }

        {

          PTSet writers({InteriorEntity,BorderEntity});
          PTSet readers({InteriorEntity,BorderEntity});

          CommunicationTestDataHandle<GV> dh(gv,codim,writers,readers,coords,periodicMap);
          gv.communicate(dh,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);
          dh.verify(
            codim,
            [](PartitionType partition, bool /* write_allowed */, bool /* read_allowed */, std::size_t writes, std::size_t reads)
            {
              if (partition == BorderEntity)
                return (writes > 0) && (reads > 0) && (writes == reads);
              else
                return (writes == 0) && (reads == 0);
            });
        }
      }

      // Need a forward declaration here
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim, const PeriodicMap<GV>& periodicMap);

      // Statically iterate over all codimensions
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV, Codim<cd>, const PeriodicMap<GV>&, std::true_type, std::true_type) {}
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV, Codim<cd>, const PeriodicMap<GV>&, std::false_type, std::true_type) {}

      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim, const PeriodicMap<GV>& periodicMap, std::true_type, std::false_type) {
        check_communication_correctness_do(gv, codim, periodicMap);
        check_communication_correctness_iter (gv, Codim<cd+1>(), periodicMap);
      }
      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim, const PeriodicMap<GV>& periodicMap, std::false_type, std::false_type) {
        check_communication_correctness_iter (gv, Codim<cd+1>(), periodicMap);
      }

      template<typename GV, int cd>
      void check_communication_correctness_iter(GV gv, Codim<cd> codim, const PeriodicMap<GV>& periodicMap) {
        check_communication_correctness_iter (gv, codim, periodicMap,
              std::integral_constant<bool, Dune::Capabilities::hasEntity<typename GV::Grid, cd>::v>(),
              std::integral_constant<bool, cd == GV::dimension + 1>());
      }

      /** \brief Check the correctness of the communication on a grid view, starting with codim 0
       *
       * \param periodicMap maps each coordinate to its representative in the periodic
       *                    unit cell, see dune/grid/test/periodicmap.hh
       */
      template<typename GV>
      void check_communication_correctness(GV gv, const PeriodicMap<GV>& periodicMap = noPeriodicMap<GV>())
      {
        check_communication_correctness_iter (gv, Codim<0>(), periodicMap);
      }

  } // namespace GridCheck
} // namespace Dune

#endif
