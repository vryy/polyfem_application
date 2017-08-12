//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 21, 2016 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "polyfem_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosPolyFEMApplication)
    {

        class_<KratosPolyFEMApplication, KratosPolyFEMApplication::Pointer,
               bases<KratosApplication>, boost::noncopyable>
               ("KratosPolyFEMApplication");

        PolyFEMApplication_AddCustomUtilitiesToPython();

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON

