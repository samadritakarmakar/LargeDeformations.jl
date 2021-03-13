module LargeDeformations
using LinearAlgebra, ForwardDiff
include("deformationTools.jl")
include("convertTools.jl")
include("hyperElasticModels.jl")


#from convertTools.jl
export getMandelIndex, convert2DTensorToMandel, convert4DTensorToMandel,
    convert2DMandelToTensor, convert4DMandelToTensor

#from deformationTools.jl
export getDeformationGradient, getDeformationGradient!,
    getJacobianDeformationGradient,
    getRightCauchyDeformation, getRightCauchyDeformation!,
    getLeftCauchyDeformation, getLeftCauchyDeformation!,
    getGreenStrain, getGreenStrain!,
    convertGreen2AlmansiStrain, convertGreen2AlmansiStrain!,
    getAlmansiStrain,
    getInvariants, getPrincipalStretches,
    convertMaterialTangent2SpatialTangent!,
    convertSecondPiola2CauchyStress, convertSecondPiola2CauchyStress!

#from hyperElasticModels.jl
export hyperElasticModel, saintVenantCauchyStress, saintVenantTangent!,
    saintVenantSpatialTangent!
end # module
