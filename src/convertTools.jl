using LinearAlgebra

function getMandelIndex(i,j)
    return 3*(j-1)+i
end

function convert2DTensorToMandel(tensor2D::Array{T, 2}) where T
    mandel2D = Array{T, 1}(undef, size(tensor2D, 1)*size(tensor2D,2))
    for j ∈ 1:size(tensor2D, 2)
        for i ∈ 1:size(tensor2D, 1)
            mandel2D[3*(j-1)+i] = tensor2D[i,j]
        end
    end
    return mandel2D
end

function convert4DTensorToMandel(tensor4D::Array{T, 4}) where T
    mandel4D = Array{T, 2}(undef, size(tensor4D, 1)*size(tensor4D,2), size(tensor4D, 3)*size(tensor4D,4))
    counter = 1
    for l ∈ 1:size(tensor4D, 4)
        for k ∈ 1:size(tensor4D, 3)
            for j ∈ 1:size(tensor4D, 2)
                for i ∈ 1:size(tensor4D, 1)
                    mandel4D[3*(j-1)+i, 3*(l-1)+k] = tensor4D[i,j,k,l]
                end
            end
        end
    end
    return mandel4D
end

function convert2DMandelToTensor(mandel2D::Array{T,1}) where T
    tensor2D = Array{T, 2}(undef, 3,3)
    for j ∈ 1:3
        for i ∈ 1:3
            tensor2D[i,j] = mandel2D[3*(j-1)+i]
        end
    end
    return tensor2D
end

function convert4DMandelToTensor(mandel4D::Array{T,2}) where T
    tensor4D = Array{T, 4}(undef, 3, 3, 3, 3)
    for l ∈ 1:3
        for k ∈ 1:3
            for j ∈ 1:3
                for i ∈ 1:3
                    tensor4D[i,j,k,l] = mandel4D[3*(j-1)+i, 3(l-1)+k]
                end
            end
        end
    end
    return tensor4D
end
