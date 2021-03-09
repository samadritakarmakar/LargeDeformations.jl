struct hyperElasticModel
    secondPiolaStress::Function
    materialTangentTensor!::Function
    cauchyStress::Function
    spatialTangentTensor!::Function
end

####Saint Venant#############

function saintVenantSecondPiolaStress(E_tensor::Array{T, 2}, λ::Float64, μ::Float64) where T
    #E_tensor = getGreenStrain(F_tensor)
    S = zeros(T, 3, 3)
    trace_E = LinearAlgebra.tr(E_tensor)
    for J ∈ 1:3
        for I ∈ 1:3
            S[I,J] += λ*δ(I,J)*trace_E + 2*μ*E_tensor[I,J]
        end
    end
    return S
end

function saintVenantSecondPiolaStress(E_mandel::Array{T, 1}, λ::Float64, μ::Float64) where T
    S = zeros(T, 9)
    #E_mandel = getGreenStrain(F_mandel)
    #E_tensor = convert2DMandelToTensor(E_mandel)

    trace_E = 0.0
    for K ∈ 1:3
        KK = getMandelIndex(K,K)
        trace_E += E_mandel[KK]
    end
    for J ∈ 1:3
        for I ∈ 1:3
            IJ = getMandelIndex(I,J)
            S[IJ] += λ*δ(I,J)*trace_E + 2*μ*E_mandel[IJ]
        end
    end
    return S
end

function saintVenantSecondPiolaStress(E::Array{T,N}, λ_μ::Tuple{Float64, Float64}) where {T, N}
    return saintVenantSecondPiolaStress(E, λ_μ[1], λ_μ[2])
end


function saintVenantCauchyStress(F::Array{T,N}, λ_μ::Tuple{Float64, Float64}) where {T, N}
    λ = λ_μ[1]
    μ = λ_μ[2]
    E = getGreenStrain(F)
    S = saintVenantSecondPiolaStress(E, λ, μ)
    return convertSecondPiola2CauchyStress(S, F)
end

function saintVenantTangent!(ℂ::Array{T, 4}, E_tensor::Array{T,2}, λ_μ::Tuple{Float64, Float64}) where T

    λ = λ_μ[1]
    μ = λ_μ[2]
    #ℂ = zeros(T, 3, 3, 3, 3)
    fill!(ℂ, 0.0)
    for L ∈ 1:3
        for K ∈ 1:3
            δ_KL = δ(K, L)
            for J ∈ 1:3
                δ_JL = δ(J,L)
                δ_JK = δ(J,K)
                for I ∈ 1:3
                    #δ_IJ = δ(I, J)
                    #δ_IK = δ(I, K)
                    #δ_IL = δ(I, L)
                    ℂ[I,J,K,L] += λ*δ(I, J)*δ_KL + μ*(δ(I, K)*δ_JL+δ_JK*δ(I, L))
                end
            end
        end
    end
    return ℂ
end

function saintVenantTangent!(ℂ::Array{T, 2}, E_mandel::Array{T,1}, λ_μ::Tuple{Float64, Float64}) where T
    λ = λ_μ[1]
    μ = λ_μ[2]
    fill!(ℂ, 0.0)
    for L ∈ 1:3
        for K ∈ 1:3
            KL = getMandelIndex(K, L)
            δ_KL = δ(K, L)
            for J ∈ 1:3
                δ_JL = δ(J,L)
                δ_JK = δ(J,K)
                for I ∈ 1:3
                    IJ = getMandelIndex(I, J)
                    #δ_IJ = δ(I, J)
                    #δ_IK = δ(I, K)
                    #δ_IL = δ(I, L)
                    ℂ[IJ, KL] += λ*δ(I, J)*δ_KL + μ*(δ(I, K)*δ_JL+δ_JK*δ(I, L))
                end
            end
        end
    end
    return ℂ
end

function saintVenantSpatialTangent!(𝕔::Array{T,4}, F::Array{T,2}, λ_μ::Tuple{Float64, Float64}) where {T, N}
    E = getGreenStrain(F)
    ℂ = zeros(T, 3,3,3,3)
    fill!(𝕔, 0.0)
    saintVenantTangent!(ℂ, E, λ_μ)
    return convertMaterialTangent2SpatialTangent!(𝕔, ℂ, F)
end

function saintVenantSpatialTangent!(𝕔::Array{T,2}, F::Array{T,1}, λ_μ::Tuple{Float64, Float64}) where {T, N}
    E = getGreenStrain(F)
    ℂ = zeros(T, 9,9)
    fill!(𝕔, 0.0)
    saintVenantTangent!(ℂ, E, λ_μ)
    return convertMaterialTangent2SpatialTangent!(𝕔, ℂ, F)
end
##### Definition of Saint Venant Hyper Elastic Model
saintVenantModel  = hyperElasticModel(saintVenantSecondPiolaStress,
    saintVenantTangent!, saintVenantCauchyStress,
    saintVenantSpatialTangent!)
