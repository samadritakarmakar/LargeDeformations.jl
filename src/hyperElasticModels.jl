struct hyperElasticModel
    cauchyStress::Function
    spatialTangentTensor::Function
end

####Saint Venant#############

function saintVenantSecondPiolaStress(E_tensor::Array{T, 2}, λ::Float64, μ::Float64) where T
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
    E_tensor = convert2DMandelToTensor(E_mandel)
    trace_E = LinearAlgebra.tr(E_tensor)
    for J ∈ 1:3
        for I ∈ 1:3
            IJ = getMandelIndex(I,J)
            S[IJ] += λ*δ(I,J)*trace_E + 2*μ*E_mandel[IJ]
        end
    end
    return S
end

function saintVenantCauchyStress(F::Array{T,N}, λ_μ::Tuple{Float64, Float64}) where {T, N}
    λ = λ_μ[1]
    μ = λ_μ[2]
    E = getGreenStrain(F)
    S = saintVenantSecondPiolaStress(E, λ, μ)
    return convertSecondPiola2CauchyStress(S, F)
end

function saintVenantTangent(F_tensor::Array{T,2}, λ_μ::Tuple{Float64, Float64}) where T
    λ = λ_μ[1]
    μ = λ_μ[2]
    ℂ = zeros(T, 3, 3, 3, 3)
    for K ∈ 1:3
        for L ∈ 1:3
            for J ∈ 1:3
                for I ∈ 1:3
                    ℂ[I,J,K,L] += λ*δ(I, J)*δ(K, L) + μ*(δ(I,K)*δ(J,L)+δ(J,K)*δ(I,L))
                end
            end
        end
    end
    return ℂ
end

function saintVenantTangent(F_mandel::Array{T,1}, λ_μ::Tuple{Float64, Float64}) where T
    λ = λ_μ[1]
    μ = λ_μ[2]
    ℂ = zeros(T, 9, 9)
    for K ∈ 1:3
        for L ∈ 1:3
            KL = getMandelIndex(K, L)
            for J ∈ 1:3
                for I ∈ 1:3
                    IJ = getMandelIndex(I, J)
                    ℂ[IJ, KL] += λ*δ(I, J)*δ(K, L) + μ*(δ(I,K)*δ(J,L)+δ(J,K)*δ(I,L))
                end
            end
        end
    end
    return ℂ
end

function saintVenantSpatialTangent(F::Array{T,N}, λ_μ::Tuple{Float64, Float64}) where {T, N}
    ℂ = saintVenantTangent(F, λ_μ)
    return convertMaterialTangent2SpatialTangent(ℂ, F)
end

##### Definition of Saint Venant Hyper Elastic Model
saintVenantModel  = hyperElasticModel(saintVenantCauchyStress, saintVenantSpatialTangent)
