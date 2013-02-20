message(AUTHOR_WARNING "CMake will only find the most current UG")

find_package(ug 3.9.1)
message(AUTHOR_WARNING "We need to test for the patch level, too!")

function(add_dune_ug_flags)
  if(UG_FOUND)
    cmake_parse_arguments(ADD_UG "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_UG_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
      set_property(DIRECTORY APPEND PROPERTY INCLUDE_DIRECTORIES ${UG_INCLUDES})
    else()
      if(ADD_UG_OBJECT)
	set(_object OBJECT)
      else(ADD_UG_OBJECT)
      set(_prefix TARGET)
	foreach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
	  target_link_libraries(${_target} ${UG_LIBRARIES} ${DUNE_LIBS})
	  set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS} APPEND
	    PROPERTY LIBRARY_FLAGS ${UG_LIBRARY_FLAGS})
	endforeach(_target ${ADD_UG_UNPARSED_ARGUMENTS})
      endif(ADD_UG_OBJECT)
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS} APPEND
	PROPERTY INCLUDE_DIRECTORIES ${UG_INCLUDES})
    endif()

    set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS} APPEND PROPERTY COMPILE_DEFINITIONS ENABLE_UGGRID ${UG_COMPILE_FLAGS})
    if(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
      set_property(${_prefix} ${ADD_UG_UNPARSED_ARGUMENTS} APPEND PROPERTY LINK_LIBRARIES ${UG_LIBRARIES} dunegrid ${DUNE_LIBS})
    endif(NOT (ADD_UG_SOURCE_ONLY OR ADD_UG_OBJECT))
  endif(UG_FOUND)
endfunction(add_dune_ug_flags)