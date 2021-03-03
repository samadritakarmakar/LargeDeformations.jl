struct hyperElasticModel
    cauchyStress::Function
    spatialTangentTensor::Function
end

####Saint Venant#############
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

function saintVenantCauchyStress(F_mandel::Array{T,1}, λ_μ::Tuple{Float64, Float64}) where T
    λ = λ_μ[1]
    μ = λ_μ[2]
    E_mandel = getGreenStrain(F_mandel)
    S_mandel = saintVenantSecondPiolaStress(E_mandel, λ, μ)
    return convertSecondPiola2CauchyStress(S_mandel, F_mandel)
end

function saintVenantSpatialTangent(F_mandel::Array{T,1}, λ_μ::Tuple{Float64, Float64}) where T
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
    return convertMaterialTangent2SpatialTangent(ℂ, F_mandel)
end

##### Definition of Saint Venant Hyper Elastic Model
saintVenantModel  = hyperElasticModel(saintVenantCauchyStress, saintVenantSpatialTangent)
