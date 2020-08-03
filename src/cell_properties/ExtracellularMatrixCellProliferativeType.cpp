#include "ExtracellularMatrixCellProliferativeType.hpp"

ExtracellularMatrixCellProliferativeType::ExtracellularMatrixCellProliferativeType()
    : AbstractCellProliferativeType(5)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ExtracellularMatrixCellProliferativeType)
