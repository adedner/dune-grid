// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_VIRTUALIZEDGRID_CAST_HH
#define DUNE_VIRTUALIZEDGRID_CAST_HH

namespace Dune {

  namespace _Utility
  {
    template <class T>
    constexpr std::string_view type_name()
    {
      using namespace std;
#ifdef __clang__
      string_view p = __PRETTY_FUNCTION__;
      return string_view(p.data() + 50, p.size() - 50 - 1);
#elif defined(__GNUC__)
      string_view p = __PRETTY_FUNCTION__;
#  if __cplusplus < 201402
      return string_view(p.data() + 52, p.size() - 52 - 1);
#  else
      return string_view(p.data() + 65, p.find(';', 65) - 65);
#  endif
#elif defined(_MSC_VER)
      string_view p = __FUNCSIG__;
      return string_view(p.data() + 100, p.size() - 100 - 7);
#endif
    }
  }  // namespace _Utility

  template<class Implementation, class Interface>
  auto&& upcast(Interface&& interface)
  {
    try
    {
      return static_cast<const Implementation&>(*interface.impl().impl_.get()).impl();
    }
    catch(const std::bad_cast& e)
    {
      std::cout << e.what() << std::endl;
      std::cout << _Utility::type_name<decltype(*interface.impl().impl_.get())>() << std::endl;
      std::cout << _Utility::type_name<const Implementation&>() << std::endl;
      DUNE_THROW(InvalidStateException, "Dynamic cast failed!");
    }
  }

}  // namespace Dune

#endif
