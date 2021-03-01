module LargeDeformations
include("deformationTools.jl")
include("convertTools.jl")
include("hyperElasticModels.jl")

#from convertTools.jl
export getMandelIndex, convert2DTensorToMandel, convert4DTensorToMandel,
    convert2DMandelToTensor, convert4DMandelToTensor

#from deformationTools.jl
export getDeformationGradient, getJacobianDeformationGradient,
    getRightCauchyDeformation, getLeftCauchyDeformation,
    getLagrangeStrain, getGreenStrain, convertGreen2AlmansiStrain,
    getAlmansiStrain, getInvariants, getPrincipalStretches,
    convertMaterialTangent2SpatialTangent, convertSecondPiola2CauchyStress

#from hyperElasticModels.jl
export hyperElasticModel, saintVenantCauchyStress,
    saintVenantSpatialTangent
end # module
