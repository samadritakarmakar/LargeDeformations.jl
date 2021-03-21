using LargeDeformations, LinearAlgebra

function hyperElasticTest()
    #∂u_∂X = zeros(9)
    #∂u_∂X[1] = 1e-4
    ∂u_∂X = rand(9)*1e-6
    F = LargeDeformations.getDeformationGradient(∂u_∂X)
    J = LargeDeformations.getJacobianDeformationGradient(F)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    #println("∂u_∂X = ", ∂u_∂X)
    #println("F_tensor = ", F_tensor)
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    λ_μ = (λ, μ)
    #println("λ_μ = ", λ_μ)
    𝔼 = LargeDeformations.getGreenStrain(F)
    #println("E = ", 𝔼)
    σ = zeros(9)
    hyperModel = LargeDeformations.saintVenantModel
    #hyperModel = LargeDeformations.neoHookeanCompressibleModel
    hyperModel.cauchyStress!(σ, F, λ_μ)
    ##############################
    S_check1 = zeros(9)
    hyperModel.secondPiolaStress!(S_check1, 𝔼, λ_μ)
    #############################
    #println("σ = ", σ)
    F_tensor = LargeDeformations.convert2DMandelToTensor(F)
    ∂u_∂X_tensor = LargeDeformations.convert2DMandelToTensor(∂u_∂X)
    ℂ = zeros(9,9)
    hyperModel.materialTangentTensor!(ℂ, 𝔼, λ_μ)
    ############################
    S_check2 = ℂ*𝔼
    ###########################
    S = LargeDeformations.convert2DMandelToTensor(S_check2)
    println("Second Piola Stress Check ", norm(S_check2- S_check1))
    σ_check = LargeDeformations.convert2DTensorToMandel(1/J*F_tensor*S*F_tensor')
    𝕔_mandel = zeros(9,9)
    hyperModel.spatialTangentTensor!(𝕔_mandel, F, λ_μ)
    #println(𝕔_mandel)
    𝕔_tensor = zeros(3,3,3,3)
    #hyperModel.spatialTangentTensor!(𝕔_tensor, F_tensor, λ_μ)
    #println(LargeDeformations.convert4DTensorToMandel(𝕔_tensor))
    println("Spatial Tangent Tensor Check ", norm(𝕔_mandel - LargeDeformations.convert4DTensorToMandel(𝕔_tensor))<1e-9)
end

function hyperElasticTestNew()
    #∂u_∂X = zeros(9)
    #∂u_∂X[1] = 1e-4
    E::Float64 = 200e3 #MPa
    ν::Float64 = 0.3
    λ = (ν*E)/((1+ν)*(1-2*ν))
    μ = E/(2*(1+ν))
    λ_μ = (λ, μ)
    ∂u_∂X_total::Array{Float64, 1} = [1e-3, zeros(Float64,8)...]
    totalSteps = 1000
    𝔼_lastStep = zeros(9)
    S_check2 = zeros(9)
    ℂ = zeros(9,9)

    for step ∈ 1:totalSteps
        ∂u_∂X = (step/totalSteps)*∂u_∂X_total
        F = LargeDeformations.getDeformationGradient(∂u_∂X)
        J = LargeDeformations.getJacobianDeformationGradient(F)
        F_tensor = LargeDeformations.convert2DMandelToTensor(F)
        #println("∂u_∂X = ", ∂u_∂X)
        #println("F_tensor = ", F_tensor)

        #println("λ_μ = ", λ_μ)
        𝔼_step = LargeDeformations.getGreenStrain(F)
        #println("E = ", 𝔼_step)
        σ = zeros(9)
        #hyperModel = LargeDeformations.saintVenantModel
        hyperModel = LargeDeformations.neoHookeanCompressibleModel
        hyperModel.cauchyStress!(σ, F, λ_μ)
        ##############################
        S_check1 = zeros(9)
        hyperModel.secondPiolaStress!(S_check1, 𝔼_step, λ_μ)
        #############################
        #println("σ = ", σ)
        F_tensor = LargeDeformations.convert2DMandelToTensor(F)
        ∂u_∂X_tensor = LargeDeformations.convert2DMandelToTensor(∂u_∂X)

        hyperModel.materialTangentTensor!(ℂ, 𝔼_step, λ_μ)
        ############################
        #println("ℂ = ", ℂ)
        S_check2 += ℂ*(𝔼_step-𝔼_lastStep)
        ###########################
        S = LargeDeformations.convert2DMandelToTensor(S_check2)
        println("Second Piola Stress Check ", norm(S_check2- S_check1))
        𝔼_lastStep = deepcopy(𝔼_step)
        #σ_check = LargeDeformations.convert2DTensorToMandel(1/J*F_tensor*S*F_tensor')
        #𝕔_mandel = zeros(9,9)
        #hyperModel.spatialTangentTensor!(𝕔_mandel, F, λ_μ)
        #println(𝕔_mandel)
        #𝕔_tensor = zeros(3,3,3,3)
        #hyperModel.spatialTangentTensor!(𝕔_tensor, F_tensor, λ_μ)
        #println(LargeDeformations.convert4DTensorToMandel(𝕔_tensor))

        #println("Spatial Tangent Tensor Check ", norm(𝕔_mandel - LargeDeformations.convert4DTensorToMandel(𝕔_tensor))<1e-9)
    end
end
