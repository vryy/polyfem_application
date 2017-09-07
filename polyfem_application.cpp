//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Jan 21, 2016$
//   Revision:            $Revision: 1.0 $
//
// 


// System includes


// External includes


// Project includes
#include "polyfem_application.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"

namespace Kratos
{
    // create the application variables here
    KRATOS_CREATE_VARIABLE(int, POLYTREE_DEBUG_LEVEL)
    KRATOS_CREATE_VARIABLE(double, MERGE_PARAMETER)

    // constructor
    KratosPolyFEMApplication::KratosPolyFEMApplication()
    {}

    // register the application to the Kratos kernel
    void KratosPolyFEMApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosPolyFEMApplication... " << std::endl;

        // register variable
        KRATOS_REGISTER_VARIABLE(POLYTREE_DEBUG_LEVEL)
        KRATOS_REGISTER_VARIABLE(MERGE_PARAMETER)
    }

} // namespace Kratos

