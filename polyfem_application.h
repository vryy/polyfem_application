//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Jan 21, 2016 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_POLYFEM_APPLICATION_H_INCLUDED)
#define KRATOS_POLYFEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/legacy_structural_app_vars.h"
#include "includes/kratos_application.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

    ///@name Kratos Globals
    ///@{

    // Variables definition
    KRATOS_DEFINE_APPLICATION_VARIABLE(POLYFEM_APPLICATION, int, POLYTREE_DEBUG_LEVEL)
    KRATOS_DEFINE_APPLICATION_VARIABLE(POLYFEM_APPLICATION, double, MERGE_PARAMETER)

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// Short class definition.
    /** Detail class definition.
    */
    class KRATOS_API(POLYFEM_APPLICATION) KratosPolyFEMApplication : public KratosApplication
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosMultiphaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(KratosPolyFEMApplication);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        KratosPolyFEMApplication();

        /// Destructor.
        ~KratosPolyFEMApplication() override {}

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void Register() override;

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            return "Application for polygonal Finite Element Method";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << Info();
            PrintData(rOStream);
        }

        ///// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            rOStream << "in KratosPolyFEMApplication:";
            KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
            rOStream << "Variables:" << std::endl;
            KratosComponents<VariableData>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Elements:" << std::endl;
            KratosComponents<Element>().PrintData(rOStream);
            rOStream << std::endl;
            rOStream << "Conditions:" << std::endl;
            KratosComponents<Condition>().PrintData(rOStream);
        }


        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{


        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{


        /// Assignment operator.
        KratosPolyFEMApplication& operator=(KratosPolyFEMApplication const& rOther);

        /// Copy constructor.
        KratosPolyFEMApplication(KratosPolyFEMApplication const& rOther);


        ///@}

    }; // Class KratosPolyFEMApplication

    ///@}


    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}


} // namespace Kratos

#endif // KRATOS_POLYFEM_APPLICATION_H_INCLUDED defined
