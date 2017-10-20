#ifndef DUNE_PYTHON_GRID_LOCALVIEW_HH
#define DUNE_PYTHON_GRID_LOCALVIEW_HH

#include <map>

#include <dune/common/visibility.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace CorePy
  {

    namespace detail
    {

      // LocalViewRegistry
      // -----------------

      template< class LocalView, class Context >
      struct LocalViewRegistry
      {
        void bind ( pybind11::handle localView, pybind11::object context )
        {
          localView.template cast< LocalView & >().bind( context.template cast< const Context & >() );
          find( localView ) = context;
        }

        void unbind ( pybind11::handle localView )
        {
          localView.template cast< LocalView & >().unbind();
          find( localView ) = pybind11::object();
        }

      private:
        pybind11::object &find ( pybind11::handle localView )
        {
          auto result = binds_.emplace( localView.ptr(), pybind11::object() );
          const auto pos = result.first;
          if( result.second )
          {
            pybind11::cpp_function remove_bind( [ this, pos ] ( pybind11::handle weakref ) {
                binds_.erase( pos );
                weakref.dec_ref();
              } );
            pybind11::weakref weakref( localView, remove_bind );
            weakref.release();
          }
          return pos->second;
        }

        std::map< void *, pybind11::object > binds_;
      };



      // localViewRegsitry
      // -----------------

      template< class LocalView, class Context >
      DUNE_EXPORT inline LocalViewRegistry< LocalView, Context > &localViewRegistry ()
      {
        // TODO add to a python module?
        static LocalViewRegistry< LocalView, Context > registry;
        return registry;
      }

    } // namespace detail



    // registerLocalView
    // -----------------

    template< class Context, class LocalView, class... options >
    void registerLocalView ( pybind11::class_< LocalView, options... > cls )
    {
      using pybind11::operator""_a;

      auto &registry = detail::localViewRegistry< LocalView, Context >();
      cls.def( "bind", [ &registry ] ( pybind11::handle self, pybind11::object context ) { registry.bind( self, context ); }, "context"_a );
      cls.def( "unbind", [ &registry ] ( pybind11::handle self ) { registry.unbind( self ); } );
    }

  } // namespace CorePy

} // namespace Dune

#endif // #ifndef DUNE_PYTHON_GRID_LOCALVIEW_HH
