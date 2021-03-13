using LargeDeformations, LinearAlgebra

function saintVenantTest()
    #∂u_∂X = zeros(9)
    #∂u_∂X[1] = 1e-4
    ∂u_∂X = rand(9)*1e-3
    F = LargeDeformations.getDeformationGradient(∂u_∂X)
    J = LargeDeformations.getJacobianDeformationGradient(F)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    println("∂u_∂X = ", ∂u_∂X)
    println("F_tensor = ", F_tensor)
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    λ_μ = (λ, μ)
    println("λ_μ = ", λ_μ)
    𝔼 = LargeDeformations.getGreenStrain(F)
    println("E = ", 𝔼)
    σ = zeros(9)
    LargeDeformations.saintVenantModel.cauchyStress!(σ, F, λ_μ)
    S_check1 = zeros(9)
    LargeDeformations.saintVenantModel.secondPiolaStress!(S_check1, 𝔼, λ_μ)
    println("σ = ", σ)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    ∂u_∂X_tensor = LargeDeformations.convert2DMandelToTensor(∂u_∂X)
    ℂ = zeros(9,9)
    LargeDeformations.saintVenantTangent!(ℂ, 𝔼, λ_μ)
    S_check2 = ℂ*𝔼
    S = LargeDeformations.convert2DMandelToTensor(S_check2)
    println("Second Piola Stress Check ", norm(S_check2-S_check1)<1e-12)
    σ_check = LargeDeformations.convert2DTensorToMandel(1/J*F_tensor*S*F_tensor')
    𝕔_mandel = zeros(9,9)
    LargeDeformations.saintVenantModel.spatialTangentTensor!(𝕔_mandel, F, λ_μ)
    #println(𝕔_mandel)
    𝕔_tensor = zeros(3,3,3,3)
    LargeDeformations.saintVenantModel.spatialTangentTensor!(𝕔_tensor, F_tensor, λ_μ)
    #println(LargeDeformations.convert4DTensorToMandel(𝕔_tensor))
    println("Spatial Tangent Tensor Check ", norm(𝕔_mandel - LargeDeformations.convert4DTensorToMandel(𝕔_tensor))<1e-9)
end
