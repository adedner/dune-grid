// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_VIRTUALIZEDGRID_TYPEERASURE_TYPEERASURE_HH
#define DUNE_GRID_VIRTUALIZEDGRID_TYPEERASURE_TYPEERASURE_HH

#include <typeinfo>

#include <dune/common/typeutilities.hh>

#include <dune/grid/virtualizedgrid/typeerasure/interfaces.hh>
#include <dune/grid/virtualizedgrid/typeerasure/polymorphicsmallobject.hh>

namespace Dune {
namespace TypeErasure {
namespace Impl {



/**
 * \brief The internal wrapper interface for type erasure
 *
 * \ingroup TypeErasure
 *
 * This class adds some foundation interfaces needed
 * for proper dynamic polymorphism and type erasure.
 *
 * The actual application interface has to be provided
 * by the user.
 *
 * \tparam Interface Class definining the internal abstract virtual interface
 */
template<class Interface>
class TypeErasureWrapperInterface :
  public Interface,
  public PolymorphicType<TypeErasureWrapperInterface<Interface>>
{
public:
  virtual const std::type_info& target_type() const = 0;
};



/**
 * \brief Base implementation of the internal wrapper interface
 *
 * \ingroup TypeErasure
 *
 * This class derives from the foundation interfaces
 * and the user defined interfaces provided by the interface
 * parameter. It will store any suitable type T to do
 * the type erasure.
 *
 * The implementation of the foundation and user interface
 * is provided by classed derived of this one.
 *
 * \tparam Interface Class definining the internal abstract virtual interface
 * \tparam T A type modelling the desired interface
 */
template<class Interface, class T>
class TypeErasureWrapperBase :
  public TypeErasureWrapperInterface<Interface>
{
public:
  using Wrapped = T;

  template<class TT, disableCopyMove<TypeErasureWrapperBase, TT> = 0>
  TypeErasureWrapperBase(TT&& t) :
    wrapped_(std::forward<TT>(t))
  {}

  //! Get mutable reference stored object
  T& get()
  {
    return wrapped_;
  }

  //! Get reference stored object
  const T& get() const
  {
    return wrapped_;
  }

protected:
  Wrapped wrapped_;
};


template<class Interface, class T>
class TypeErasureWrapperBase<Interface, std::reference_wrapper<T>> :
  public TypeErasureWrapperInterface<Interface>
{
public:
  using Wrapped = T;

  template<class TT, disableCopyMove<TypeErasureWrapperBase, TT> = 0>
  TypeErasureWrapperBase(TT&& t) :
    wrapped_(std::forward<TT>(t))
  {}

  //! Get mutable reference stored object
  T& get()
  {
    return wrapped_.get();
  }

  //! Get reference stored object
  const T& get() const
  {
    return wrapped_.get();
  }

protected:
  std::reference_wrapper<T> wrapped_;
};


/**
 * \brief Implementation of the internal wrapper interface
 *
 * \ingroup TypeErasure
 *
 * This class implements the foundation and user interfaces
 * of the internal type erasure wrapper.
 *
 * The foundation interface of TypeErasureWrapperInterface is directly
 * implemented here whereas the user interface is implemented
 * by deriving from the user-provides Implementation template.
 *
 * The Implementation is a template taking one class template
 * parameter. It should directly or indirectly derive from this
 * class and inherit its constructors.
 * In order to forward the implemented methods to the erased
 * type it can use the wrapper_ member of this base class being
 * of this type.
 *
 * \tparam Interface Class defining the internal abstract virtual interface
 * \tparam Implementation Class defining the implementation of the abstract methods of Interface
 * \tparam T A type modeling the desired interface
 */
template<class Interface, template<class> class Implementation, class T>
class TypeErasureWrapperImplementation :
  public Implementation<TypeErasureWrapperBase<Interface, T> >
{
public:

  //! Construct wrapper from object
  template<class TT, disableCopyMove<TypeErasureWrapperImplementation, T> = 0>
  TypeErasureWrapperImplementation(TT&& t) :
    Implementation<TypeErasureWrapperBase<Interface, T> >(std::forward<TT>(t))
  {}

  //! Implementation of PolymorphicType::clone()
  virtual TypeErasureWrapperImplementation* clone() const
  {
    return new TypeErasureWrapperImplementation(*this);
  }

  //! Implementation of PolymorphicType::clone(void* buffer)
  virtual TypeErasureWrapperImplementation* clone(void* buffer) const
  {
    return new (buffer) TypeErasureWrapperImplementation(*this);
  }

  //! Implementation of PolymorphicType::move(void* buffer)
  virtual TypeErasureWrapperImplementation* move(void* buffer)
  {
    return new (buffer) TypeErasureWrapperImplementation(std::move(*this));
  }

  //! Get type of stored object
  virtual const std::type_info& target_type() const
  {
    return typeid(T);
  }
};

} // namespace Dune::Functions::Imp



/**
 * \brief Base class for type-erased interface wrapper
 *
 * \ingroup TypeErasure
 *
 * This is meant as a base class for the type-erased interface
 * wrapper that is actually visible to the user. By deriving
 * from this you get small object optimization for the internal
 * polymorphic wrapper.
 */
template<class Interface, template<class> class Implementation, size_t bufferSize = 56>
class TypeErasureBase
{
public:
  template <class T>
  using WrapperImplementation = Impl::TypeErasureWrapperImplementation<Interface, Implementation, T>;

  //! Construct wrapper from object
  template<class T, disableCopyMove<TypeErasureBase, T> = 0 >
  TypeErasureBase(T&& t) :
    wrapped_(WrapperImplementation<std::decay_t<T>>(std::forward<T>(t)))
  {}

  //! Default constructor
  TypeErasureBase() = default;

  //! Get mutable reference to wrapped object
  Interface& asInterface()
  {
    return wrapped_.get();
  }

  //! Get reference to wrapped object
  const Interface& asInterface() const
  {
    return wrapped_.get();
  }

  //! Get type of stored object
  const std::type_info& target_type() const
  {
    return wrapped_.get().target_type();
  }

protected:
  PolymorphicSmallObject<Impl::TypeErasureWrapperInterface<Interface>, bufferSize > wrapped_;
};


}} // namespace Dune::Typeerasure



#endif // DUNE_GRID_VIRTUALIZEDGRID_TYPEERASURE_TYPEERASURE_HH
