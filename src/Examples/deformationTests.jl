using LargeDeformations

function saintVenantTest()
    ∂u_∂X = zeros(9)
    ∂u_∂X[1] = 1e-4
    F = LargeDeformations.getDeformationGradient(∂u_∂X)
    J = LargeDeformations.getJacobianDeformationGradient(F)
    LargeDeformations.convert2DMandelToTensor(F)
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    λ_μ = (λ, μ)
    println("λ_μ = ", λ_μ)
    𝔼 = LargeDeformations.getGreenStrain(F)
    println("E = ", 𝔼)
    σ = LargeDeformations.saintVenantModel.cauchyStress(F, λ_μ)
    println("σ = ", σ)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    ∂u_∂X_tensor = LargeDeformations.convert2DMandelToTensor(∂u_∂X)
    ∂u_∂x = LargeDeformations.convert2DTensorToMandel(∂u_∂X_tensor*inv(F_tensor))
    println("∂u_∂x = ", ∂u_∂x)
    𝕔 = LargeDeformations.saintVenantModel.spatialTangentTensor(F, λ_μ)
    println(𝕔)
    𝕔*inv(F*F')*E
end
