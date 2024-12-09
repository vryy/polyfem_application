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

// Project includes
#include "includes/define_python.h"
#include "polyfem_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

    BOOST_PYTHON_MODULE(KratosPolyFEMApplication)
    {

        using namespace boost::python;

        class_<KratosPolyFEMApplication, KratosPolyFEMApplication::Pointer,
               bases<KratosApplication>, boost::noncopyable>
               ("KratosPolyFEMApplication");

        PolyFEMApplication_AddCustomUtilitiesToPython();

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( POLYTREE_DEBUG_LEVEL )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( MERGE_PARAMETER )

    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON
