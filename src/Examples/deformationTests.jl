using LargeDeformations, LinearAlgebra

function saintVenantTest()
    ∂u_∂X = zeros(9)
    ∂u_∂X[1] = 1e-4
    F = LargeDeformations.getDeformationGradient(∂u_∂X)
    J = LargeDeformations.getJacobianDeformationGradient(F)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    println("F_tensor = ", F_tensor)
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    λ_μ = (λ, μ)
    println("λ_μ = ", λ_μ)
    𝔼 = LargeDeformations.getGreenStrain(F)
    println("E = ", 𝔼)
    σ = LargeDeformations.saintVenantModel.cauchyStress(F, λ_μ)
    S_check1 = LargeDeformations.saintVenantModel.secondPiolaStress(F, λ_μ)
    println("σ = ", σ)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    ∂u_∂X_tensor = LargeDeformations.convert2DMandelToTensor(∂u_∂X)
    ℂ = LargeDeformations.saintVenantTangent(F, λ_μ)
    S_check2 = ℂ*𝔼
    S = LargeDeformations.convert2DMandelToTensor(S_check2)
    println("Second Piola Stress Check ", norm(S_check2-S_check1)<1e-12)
    σ_check = LargeDeformations.convert2DTensorToMandel(1/J*F_tensor*S*F_tensor')
    𝕔_mandel = LargeDeformations.saintVenantModel.spatialTangentTensor(F, λ_μ)
    𝕔_tensor = LargeDeformations.saintVenantModel.spatialTangentTensor(F_tensor, λ_μ)
    println("Spatial Tangent Tensor Check ", 𝕔_mandel == LargeDeformations.convert4DTensorToMandel(𝕔_tensor))
end
