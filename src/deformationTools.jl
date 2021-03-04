function δ(i::Int64,j::Int64)
    return i==j ? 1.0 : 0.0
end

function getDeformationGradient(∂u_∂X_tensor::Array{T, 2}) where T
    F = zeros(T, 3,3)
    for J ∈ 1:3
        for i ∈ 1:3
            F[i,J] += δ(i,J) + ∂u_∂X_tensor[i,J]
        end
    end
    return F
end

function getDeformationGradient(∂u_∂X_mandel::Array{T, 1}) where T
    F = zeros(T, 9)
    for J ∈ 1:3
        for i ∈ 1:3
            iJ = getMandelIndex(i,J)
            F[iJ] += δ(i,J) + ∂u_∂X_mandel[iJ]
        end
    end
    return F
end

function getJacobianDeformationGradient(F_tensor::Array{T, 2}) where T
    return LinearAlgebra.det(F_tensor)
end

function getJacobianDeformationGradient(F_mandel::Array{T, 1}) where T
    F_tensor = convert2DMandelToTensor(F_mandel)
    return getJacobianDeformationGradient(F_tensor)
end

function getRightCauchyDeformation(F_tensor::Array{T, 2}) where T
    return F'*F
end

function getRightCauchyDeformation(F_mandel::Array{T, 1}) where T
    C = zeros(T, 9)
    for J ∈ 1:3
        for I ∈ 1:3
            IJ = getMandelIndex(I,J)
            for i ∈ 1:3
                iI = getMandelIndex(i,I)
                iJ = getMandelIndex(i,J)
                C[IJ] += F_mandel[iI]*F_mandel[iJ]
            end
        end
    end
    return C
end

function getLeftCauchyDeformation(F_tensor::Array{T, 2}) where T
    return F*F'
end

function getLeftCauchyDeformation(F_mandel::Array{T, 1}) where T
    b = zeros(T, 9)
    for j ∈ 1:3
        for i ∈ 1:3
            ij = getMandelIndex(i,j)
            for I ∈ 1:3
                iI = getMandelIndex(i,I)
                jI = getMandelIndex(j,I)
                b[ij] += F_madel[iI]*F_madel[jI]
            end
        end
    end
    return b
end

function getGreenStrain(F_tensor::Array{T, 2}) where T
    C = getRightCauchyDeformation(F_tensor)
    E = zeros(T, 3,3)
    for j ∈ 1:3
        for i ∈ 1:3
            E[i,j] += 0.5*(C[i,j]- δ(i,j))
        end
    end
    return E
end

function getGreenStrain(F_mandel::Array{T, 1}) where T
    C = getRightCauchyDeformation(F_mandel)
    E = zeros(T, 9)
    for j ∈ 1:3
        for i ∈ 1:3
            ij = getMandelIndex(i,j)
            E[ij] += 0.5*(C[ij]- δ(i,j))
        end
    end
    return E
end

function convertGreen2AlmansiStrain(E_tensor::Array{T, 2}, F_tensor::Array{T, 2}) where T

    Finv = inv(F_tensor)
    return Finv'*E_tensor*Finv
end

function convertGreen2AlmansiStrain(E_mandel::Array{T, 1}, F_mandel::Array{T, 1}) where T
    E_tensor = convert2DMandelToTensor(E_mandel)
    F_tensor = convert2DMandelToTensor(F_mandel)
    e_tensor = convertGreen2AlmansiStrain(E_tensor, F_tensor)
    return convert2DTensorToMandel(e_tensor)
end

function getAlmansiStrain(F_tensor::Array{T, 2}) where T
    #e = F⁻ᵀ.E.F⁻¹
    E_tensor = getGreenStrain(F_tensor)
    return convertGreen2AlmansiStrain(E_tensor, F_tensor)
end

function getAlmansiStrain(F_mandel::Array{T, 1}) where T
    F_tensor = convert2DMandelToTensor(F_mandel)
    e_tensor = getAlmansiStrain(F_tensor)
    return convert2DTensorToMandel(e_tensor)
end

function getInvariants(F_tensor::Array{T, 2}) where T
    I = LinearAlgebra.tr(F_tensor)
    II = 0.5*(I^2 - LinearAlgebra.tr(F_tensor*F_tensor))
    III = det(F_tensor)
    return I, II, III
end

function getInvariants(F_mandel::Array{T, 1}) where T
    F_tensor = convert2DMandelToTensor(F_mandel)
    return getInvariants(F_tensor)
end

function getPrincipalStretches(F_tensor::Array{T, 2}) where T
    C = getRightCauchyDeformation(F_tensor)
    λ² = LinearAlgebra.eigvals(C)
    return sqrt.(λ²)
end

function getPrincipalStretches(F_mandel::Array{T, 1}) where T
    F_tensor = convert2DMandelToTensor(F_mandel)
    return getPrincipalStretches(F_tensor)
end

function convertMaterialTangent2SpatialTangent(ℂ_tensor::Array{T1, 4}, F_tensor::Array{T2, 2}) where {T1, T2}
    J = getJacobianDeformationGradient(F_tensor)
    𝕔 = zeros(T1, 3,3,3,3)
    for l ∈ 1:3
        for k ∈ 1:3
            for j ∈ 1:3
                for i ∈ 1:3
                    for Q ∈ 1:3
                        for P ∈ 1:3
                            for N ∈ 1:3
                                for M ∈ 1:3
                                    𝕔[i,j,k,l] += F_tensor[i,M]*F_tensor[j,N]*ℂ_tensor[M,N,P,Q]*F_tensor[k,P]*F_tensor[l,Q]/J
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return 𝕔
end

function convertMaterialTangent2SpatialTangent(ℂ::Array{T1, 2}, F_mandel::Array{T2, 1}) where {T1, T2}
    J = getJacobianDeformationGradient(F_mandel)
    𝕔 = zeros(T1, 9, 9)
    for l ∈ 1:3
        for k ∈ 1:3
            kl = getMandelIndex(k,l)
            for j ∈ 1:3
                for i ∈ 1:3
                    ij = getMandelIndex(i,j)
                    for Q ∈ 1:3
                        lQ = getMandelIndex(l, Q)
                        for P ∈ 1:3
                            PQ = getMandelIndex(P, Q)
                            kP = getMandelIndex(k, P)
                            for N ∈ 1:3
                                jN = getMandelIndex(j,N)
                                for M ∈ 1:3
                                    iM = getMandelIndex(i,M)
                                    MN = getMandelIndex(M, N)
                                    𝕔[ij, kl] += F_mandel[iM]*F_mandel[jN]*ℂ[MN, PQ]*F_mandel[kP]*F_mandel[lQ]/J
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return 𝕔
end

function convertSecondPiola2CauchyStress(S_tensor::Array{T1,2}, F_tensor::Array{T2,2}) where {T1, T2}
    J = getJacobianDeformationGradient(F_tensor)
    return 1/J*(F_tensor*S_tensor*F_tensor')
end

function convertSecondPiola2CauchyStress(S_mandel::Array{T1,1}, F_mandel::Array{T2,1}) where {T1, T2}
    F_tensor = convert2DMandelToTensor(F_mandel)
    S_tensor = convert2DMandelToTensor(S_mandel)
    σ_tensor = convertSecondPiola2CauchyStress(S_tensor, F_tensor)
    return convert2DTensorToMandel(σ_tensor)
end
