/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "FibroblastStateDependentCollagenOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

FibroblastStateDependentCollagenOdeSystem::FibroblastStateDependentCollagenOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<FibroblastStateDependentCollagenOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - C (Collagen)
     * 
     */

//    SetDefaultInitialCondition(0, 0.01); // soon overwritten
//    SetDefaultInitialCondition(1, 0.01); // soon overwritten
//    SetDefaultInitialCondition(2, 0.01); // soon overwritten

    this->mParameters.push_back(0.0); // EPF status
    this->mParameters.push_back(1.0); // Production rate of collagen
    this->mParameters.push_back(0.1); // Degradation rate of collagen

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

FibroblastStateDependentCollagenOdeSystem::~FibroblastStateDependentCollagenOdeSystem()
{
}

void FibroblastStateDependentCollagenOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    // state values
    double C = rY[0]; // collagen

    // Get the parameters
    double is_epf = this->mParameters[0]; // Check whether or not cell is an EPF fibroblast
    double p_c = this->mParameters[1]; // Production rate of collagen
    double d_c = this->mParameters[2]; // Degradation rate of collagen

    // calculations
    double reaction_1 = is_epf * p_c * C;
    double reaction_2 = d_c * C;

    // ODEs
    rDY[0] = reaction_1 - reaction_2; // dC/dt = p_C*delta(EPF)*C - d_C*C
}

template<>
void CellwiseOdeSystemInformation<FibroblastStateDependentCollagenOdeSystem>::Initialise()
{

    // State variable, which is collagen amount
    this->mVariableNames.push_back("collagen");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01); // will be filled in later

    // Multiplier to check whether or not cells are EPF fibroblasts
    this->mParameterNames.push_back("epf");
    this->mParameterUnits.push_back("non-dim");

    // Production rate of collagen
    this->mParameterNames.push_back("p_c");
    this->mParameterUnits.push_back("non-dim");

    // Degradation rate of collagen
    this->mParameterNames.push_back("d_c");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FibroblastStateDependentCollagenOdeSystem)

