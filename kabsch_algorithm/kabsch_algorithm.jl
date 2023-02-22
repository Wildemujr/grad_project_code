using BioStructures
using LinearAlgebra
using Statistics
using CoordinateTransformations, Rotations, StaticArrays
using Random













# function prune_q_points(input_matrix, points_q)
    
#     pruned_matrix = Matrix{Float64}(undef, length(points_q), 3)

#     for (i,j) in enumerate(points_q)

#         # * Find the rows where the second column is equal to i
#         rows = findall(x -> x == j, input_matrix[:, 2])

#         # * Extract the submatrix with the matching rows, and then find
#         # * the minimum value in the third column of the submatrix
#         submatrix = input_matrix[rows, :]
#         min_value, min_idx = findmin(submatrix[:, 3])

#         # * Extract the values from the first and second columns of the submatrix
#         first_val = submatrix[min_idx, 1]
#         second_val = submatrix[min_idx, 2]

#         # * Add the values to the pruned matrix
#         pruned_matrix[i, :] = [first_val, second_val, min_value]

#     end

#     return pruned_matrix

# end



# function compute_cross_covariance(P::Matrix{Float64}, Q::Matrix{Float64}, gaussian_kernel::Function, bandwidth_param::Float64)

# function compute_cross_covariance(P::Matrix{Float64}, Q::Matrix{Float64})


#     # * Initialize the cross-covariance matrix
#     cross_covariance_matrix::Matrix{Float64} = zeros(Float64, size(P)[2], size(Q)[2])

#     # * Computing the geometric centroids for the two point clouds.
#     P_mean = mean(P, dims=1)
#     Q_mean = mean(Q, dims=1)
    
#     n::Int64 = max(size(P)[1], size(Q)[1])
#     m::Int64 = max(size(P)[2], size(Q)[2])
    
#     for i = 1:m
#         for j = 1:m
#             P_point = P[:, i:i]
#             Q_point = Q[:, j:j]
#             covariance = ( (P_point - P_mean) * transpose( (Q_point - Q_mean) ) ) / n |> sum           
#             cross_covariance_matrix[i,j] = covariance
#         end
#     end

#     return cross_covariance_matrix

# end


function compute_cross_covariance(P::Matrix{Float64}, Q::Matrix{Float64})


    # * Initialize the cross-covariance matrix
    cross_covariance_matrix::Matrix{Float64} = zeros(Float64, size(P)[2], size(Q)[2])




    # # * Computing the geometric centroids for the two point clouds.
    # P_mean = mean(P, dims=1)
    # Q_mean = mean(Q, dims=1)
    
    # n::Int64 = max(size(P)[1], size(Q)[1])
    # m::Int64 = max(size(P)[2], size(Q)[2])
    
    # for i = 1:m
    #     for j = 1:m
    #         P_point = P[:, i:i]
    #         Q_point = Q[:, j:j]
    #         covariance = ( (P_point - P_mean) * transpose( (Q_point - Q_mean) ) ) / n |> sum           
    #         cross_covariance_matrix[i,j] = covariance
    #     end
    # end

    return cross_covariance_matrix

end



function apply_to_structure(structure, new_coords)
    idx = 1
    for mod in structure
        for ch in mod
            for res in ch
                for at in res
                    coords!(at, new_coords[idx, :])
                    # println(at.coords)
                    idx == size(new_coords)[1] ? break : nothing
                    idx += 1
                end
            end
        end
    end
    return 0 
end


# function rodrigues_rotation(u, θ) 
#     # * https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
#     I3 = Matrix{Float64}(I, 3, 3)
#     u_cross = [0 -u[3] u[2]; u[3] 0 -u[1]; -u[2] u[1] 0]
#     R = I3 + sin(θ) * u_cross + (1 - cos(θ)) * u_cross^2
#     return R
# end


# function random_rotation()
#     # * Generate a random rotation matrix
#     θ, u = rand() * 2π, randn(3)
#     u = u ./ norm(u)
#     R = rodrigues_rotation(u, θ)
#     return R 
# end


# Function to extract coordinates from PDB file
function get_coords(pdbfile::String)
    # * Read in the structure and collect all of the atoms
    structure = read(pdbfile, PDB)
    atoms = collectatoms(structure)
    # coordinates = Array{Float64}(undef, length(atoms), 3)

    coordinates = Matrix{Float64}(undef, length(atoms), 3)

    # * Now, obtain the coordinates for each atom
    for (i, atom) in enumerate(atoms)
        coordinates[i, :] = atom.coords
    end
    
    return coordinates
end



function kabsch_alignment(pdb_coords_A, pdb_coords_B, iterations::Int64, input_structure_1, input_structure_2)

    # * Initialize the input matrix variables representing each point cloud.
    structure_A = input_structure_1
    structure_B = input_structure_2
    matrix_A = pdb_coords_A
    matrix_B = pdb_coords_B

    # * truncate the matrices to the same size
    trunc_matrix_A = matrix_A[1:min(size(matrix_A)[1], size(matrix_B)[1]), :]
    trunc_matrix_B = matrix_B[1:min(size(matrix_A)[1], size(matrix_B)[1]), :]
    
    # * Translate and center the two point clouds
    centered_matrix_A = trunc_matrix_A .- mean(trunc_matrix_A, dims=1)
    centered_matrix_B = trunc_matrix_B .- mean(trunc_matrix_B, dims=1)
    


    # total_len = size(trunc_matrix_A)[1]
    # total_idxs = collect(1:total_len)
    
    # trunc_matrix_A = [total_idxs trunc_matrix_A]
    # trunc_matrix_B = [total_idxs trunc_matrix_B]



    # # * Find centroids of the points
    # # centroid_A = mean(matrix_A, dims=1)
    # # centroid_B = mean(matrix_B, dims=1)
    # # trunc_centroid_A = mean(trunc_matrix_A, dims=1)
    # # trunc_centroid_B = mean(trunc_matrix_B, dims=1)

    # * Find variance for scaling factor
    # var_A = var(matrix_A, dims=1)
    # var_B = var(matrix_B, dims=1)
    # trunc_var_A = var(trunc_matrix_A, dims=1)
    # trunc_var_B = var(trunc_matrix_B, dims=1)

    # # * Center the two point clouds
    # centered_matrix_A = trunc_matrix_A .- trunc_centroid_A
    # centered_matrix_B = trunc_matrix_B .- trunc_centroid_B
    

    

    # best_idx = collect(1:total_len)


    # * Initialize the RMSD variables
    threshold_rmsd = 2.0
    iterations_so_far = 0
    previous_RMSD = 1e6
    best_RMSD = 1e6


    while iterations_so_far < iterations

        # * Increment the iteration counter
        if iterations_so_far % 50 == 0
            println("Iteration: $(iterations_so_far)  ;  Lowest RMSD: $(best_RMSD)") ; print("\n")
        end
        
        # * calculate the covariance matrix
        K = transpose(centered_matrix_B) * centered_matrix_A
        
        # * With the cross-correlation matrix computed, we can now compute SVD of the matrix.
        F = svd(K)

        # * Calculate the rotation matrix from the previous SVD computation.
        R = F.V * transpose(F.U)
        d = det(F.V * transpose(F.U))
        I3 = Matrix{Float64}(I, 3, 3)
        
        # * Use ternary operator to see if the determinate is negative.
        # * If so, then correct our rotation matrix to ensure a right-handed coordinate system.
        # * If the determinant is positive, then we are good to go and do nothing.
        d < 0 ? (I3[3, :] = -1 * I3[3, :]) :  nothing
        
        # * Calculate the optimal rotation matrix R, and apply to the point cloud.
        # opt_R = F.Vt * (I3 * transpose(F.U))
        opt_R = F.U * (I3 * F.Vt)
        # A_aligned = (centered_matrix_B * transpose(opt_R) ) .+ mean(centered_matrix_B, dims=1) 
        
        superimposed = (opt_R * transpose(centered_matrix_B)) .+ mean(centered_matrix_A, dims=1)
        
        rmsd = sqrt( sum( (A_aligned .- centered_matrix_B).^2 ) / size(A_aligned)[1] )
        display(rmsd)
        


        apply_to_structure(structure_B, superimposed)
        writepdb("new_P_$(iterations_so_far).pdb", structure_B)

        # apply_to_structure(structure_B, centered_matrix_B)
        # writepdb("new_Q_$(iterations_so_far).pdb", structure_B)

        # iterations_so_far += 1

    #     if rmsd < previous_RMSD && rmsd < best_RMSD

    #         best_RMSD = rmsd
    #         previous_RMSD = rmsd
            
    #         centered_matrix_A = A_aligned
    #         apply_to_structure(structure_A, A_aligned)
    #         writepdb("new_P_$(iterations_so_far).pdb", structure_A)

    #         apply_to_structure(structure_B, centered_matrix_B)
    #         writepdb("new_Q_$(iterations_so_far).pdb", structure_B)

    #         iterations_so_far = iterations_so_far + 1

    #     else
    #         previous_RMSD = rmsd
    #         centered_matrix_A = A_aligned
    #         iterations_so_far = iterations_so_far + 1
    #     end

    # end






        # if rmsd < best_RMSD
        #     best_RMSD = rmsd
        #     # best_idx = idx

        #     # * Compute the distances between equivalent pairs of atoms.
        #     d = sqrt.(sum((centered_matrix_A .- centered_matrix_B).^2, dims=2))
        #     idx = [i[1] for i in findall(d .> 8.0)] # * Find those atoms with the largest distance values.
        #     trunc_matrix_A = trunc_matrix_A[setdiff(1:end, idx), :]
        #     trunc_matrix_B = trunc_matrix_B[setdiff(1:end, idx), :]    
        #     iterations_so_far = iterations_so_far + 1

        #     # apply_to_structure(structure_A, centered_matrix_A)
        #     # writepdb("new_P_$(iterations_so_far).pdb", structure_A)
        # else
        #     # * Compute the distances between equivalent pairs of atoms.
        #     d = sqrt.(sum((centered_matrix_A .- centered_matrix_B).^2, dims=2))
        #     idx = [i[1] for i in findall(d .> 8.0)] # * Find those atoms with the largest distance values.
        #     trunc_matrix_A = trunc_matrix_A[setdiff(1:end, idx), :]
        #     trunc_matrix_B = trunc_matrix_B[setdiff(1:end, idx), :]    
        #     iterations_so_far = iterations_so_far + 1
        # end

        # if rmsd <= threshold_rmsd
        #     apply_to_structure(structure_A, centered_matrix_A)
        #     writepdb("new_P_$(iterations_so_far).pdb", structure_A)
        #     break
        # end

        # display(rmsd)
        




        # display(idx)


        # if rmsd < best_RMSD
        #     best_RMSD = rmsd
        #     best_idx = idx


        # # # * Calculate the RMSD between the two point clouds after the rotation has
        # # # * been applied.
        # rmsd = sqrt( sum( (A_aligned .- centered_matrix_B).^2 ) / size(A_aligned)[1] )
        # # display(rmsd)



    # apply_to_structure(structure_B, centered_matrix_B)
    # writepdb("new_Q_$(iterations_so_far).pdb", structure_B)

end




# function kabsch_alignment(pdb_coords_A, pdb_coords_B, iterations::Int64, input_structure_1, input_structure_2)

#     # * Initialize the input matrix variables representing each point cloud.
#     structure_A = input_structure_1
#     structure_B = input_structure_2
#     matrix_A = pdb_coords_A
#     matrix_B = pdb_coords_B

#     # * truncate the matrices to the same size
#     trunc_matrix_A = matrix_A[1:min(size(matrix_A)[1], size(matrix_B)[1]), :]
#     trunc_matrix_B = matrix_B[1:min(size(matrix_A)[1], size(matrix_B)[1]), :]
    
#     # * Find centroids of the points
#     # centroid_A = mean(matrix_A, dims=1)
#     # centroid_B = mean(matrix_B, dims=1)
#     trunc_centroid_A = mean(trunc_matrix_A, dims=1)
#     trunc_centroid_B = mean(trunc_matrix_B, dims=1)

#     # * Find variance for scaling factor
#     # var_A = var(matrix_A, dims=1)
#     # var_B = var(matrix_B, dims=1)
#     trunc_var_A = var(trunc_matrix_A, dims=1)
#     trunc_var_B = var(trunc_matrix_B, dims=1)

#     # * Center the two point clouds
#     centered_matrix_A = trunc_matrix_A .- trunc_centroid_A
#     centered_matrix_B = trunc_matrix_B .- trunc_centroid_B

#     iterations_so_far = 0
#     previous_RMSD = 1e6
#     best_RMSD = 1e6

#     while iterations_so_far < iterations
        
#         # # * Increment the iteration counter
#         if iterations_so_far % 100 == 0
#             println("Iteration: $(iterations_so_far)  ;  Lowest RMSD: $(best_RMSD)") ; print("\n")
#         end


#         # # * calculate the covariance matrix
#         K = transpose(centered_matrix_A) * centered_matrix_B

#         # * With the cross-correlation matrix computed, we can now compute SVD of the matrix.
#         F = svd(K)

#         # * Calculate the rotation matrix from the previous SVD computation.
#         R = F.U * F.Vt
#         d = det(F.V * transpose(F.U))
#         I3 = Matrix{Float64}(I, 3, 3)
        
#         # * Use ternary operator to see if the determinate is negative.
#         # * If so, then correct our rotation matrix to ensure a right-handed coordinate system.
#         # * If the determinant is positive, then we are good to go and do nothing.
#         d < 0 ? (I3[3, :] = -1 * I3[3, :]) :  nothing

#         # * Calculate the optimal rotation matrix R, and apply to the point cloud.
#         # # opt_R = F.Vt * (I3 * transpose(F.U))
#         opt_R = F.V * (I3 * transpose(F.U))
#         A_aligned = centered_matrix_A * opt_R

#         # * Calculate the RMSD between the two point clouds after the rotation has
#         # * been applied.
#         rmsd = sqrt( sum( (A_aligned .- centered_matrix_B).^2 ) / size(A_aligned)[1] )
#         # display(rmsd)

#         if rmsd < previous_RMSD && rmsd < best_RMSD

#             best_RMSD = rmsd
#             previous_RMSD = rmsd
            
#             centered_matrix_A = A_aligned
#             apply_to_structure(structure_A, A_aligned)
#             writepdb("new_P_$(iterations_so_far).pdb", structure_A)

#             # apply_to_structure(structure_B, centered_matrix_B)
#             # writepdb("new_Q_$(iterations_so_far).pdb", structure_B)

#             iterations_so_far = iterations_so_far + 1

#         else
#             previous_RMSD = rmsd
#             centered_matrix_A = A_aligned
#             iterations_so_far = iterations_so_far + 1
#         end

#     end

#     apply_to_structure(structure_B, centered_matrix_B)
#     writepdb("new_Q_$(iterations_so_far).pdb", structure_B)










    #     # * See if the RMSD is less than the previous RMSD
    #     if rmsd < previous_RMSD 
    #         if rmsd < best_RMSD
    #             best_RMSD = rmsd
    #             centered_matrix_A = A_aligned

    #             apply_to_structure(structure_A, A_aligned)
    #             apply_to_structure(structure_B, centered_matrix_B)
    #             writepdb("new_P_$(iterations_so_far).pdb", structure_A)
    #             writepdb("new_Q_$(iterations_so_far).pdb", structure_B)

    #             iterations_so_far = iterations_so_far + 1

    #         else
                
    #             apply_to_structure(structure_A, A_aligned)
    #             apply_to_structure(structure_B, centered_matrix_B)
    #             writepdb("new_P_$(iterations_so_far).pdb", structure_A)
    #             writepdb("new_Q_$(iterations_so_far).pdb", structure_B)
    #             iterations_so_far = iterations_so_far + 1

    #         end

            
    #         # break
    #     else
    #         # * If not, then we need to update the previous RMSD and continue
    #         previous_RMSD = rmsd
    #         centered_matrix_A = A_aligned
    #         iterations_so_far = iterations_so_far + 1
    #     end

    #     return 0

    # end




    # apply_to_structure(structure_A, P_rotatated)
    # writepdb("new_P.pdb", structure_A)

    # apply_to_structure(structure_B, centered_matrix_B)
    # writepdb("new_Q.pdb", structure_B)


    # * Calculate the RMSD between the two point clouds after the rotation has
    # * been applied.



    # # * Initialize variables responsible for termination criteria.
    # ϵ = 1e-6
    # error_min = 1e6
    # max_iter = iterations

    # # * Initialize the iteration counter
    # iteration = 0
    # iterations_since_last_change = 0

    # for i = 1:max_iter

    #     iteration += 1
        
    #     # # * Increment the iteration counter
    #     # if i % 10 == 0
    #     #     println("Iteration: $(iteration)  ;  Error: $(error_min) ; Damping Factor: $(damping_factor)") ; print("\n\n")
    #     # end

    #     # * Increment the iteration counter
    #     if i % 10 == 0
    #         println("Iteration: $(iteration)  ;  Error: $(error_min) ; Iterations Since Last Change: $(iterations_since_last_change) ") ; print("\n\n")
    #     end

    #     # * Compute the point correspondences between the two point clouds A and B.
    #     min_dist_pts, σ = pt_correspondence(matrix_A, matrix_B)
    #     unique_source_points = unique(min_dist_pts[:, 2])
    #     point_correspondence = prune_q_points(min_dist_pts, unique_source_points)

    #     # * Here, we're isolating the indexes of the points that correspond to each other in the two point clouds.
    #     # * We then use these indexes to extract the coordinates of the points that correspond to each other.
    #     # * This will be used as input for the cross-covariance matrix computation.
    #     P_idxs::Vector{Int64} = map( x -> trunc(Int64, x), point_correspondence[:, 1])
    #     Q_idxs::Vector{Int64} = map( x -> trunc(Int64, x), point_correspondence[:, 2])

    #     # * Extracting the coordinates of the points that correspond to each other.
    #     extracted_coords_A::Matrix{Float64} = matrix_A[P_idxs, :]
    #     extracted_coords_B::Matrix{Float64} = matrix_B[Q_idxs, :]

    #     # * Feeding the coordinates into the cross-covariance matrix computation function.
    #     # * Further, isolating the geometric centroids for each point cloud.
    #     # K, P_center, Q_center = compute_cross_covariance(extracted_coords_A, extracted_coords_B)
    #     K, P_center, Q_center = compute_cross_covariance(extracted_coords_A, extracted_coords_B, gaussian_kernel,  σ)

    #     # * With the cross-correlation matrix computed, we can now compute SVD of the matrix.
    #     F = svd(K)

    #     # * Calculate the rotation matrix and translation vector from previous SVD computation.
    #     R = F.U * F.Vt
    #     t = Q_center - R * P_center

    #     # Rn = R * (1 - damping_factor) + Rn * damping_factor
    #     # tn = t * (1 - damping_factor) + tn * damping_factor
        
    #     # * Apply the rotation and translation to the point cloud A.
    #     new_P = transpose(R * transpose(matrix_A) .+ t)

    #     # * Calculate the error based on the expression for t.
    #     error = norm( Q_center - R * P_center )
    #     matrix_A = new_P

    #     # * Check if the error is less than the previous iterations.
    #     if error < error_min

    #         error_min = error
    #         iterations_since_last_change = 0

    #         # * Write out the newly transformed PDB file.
    #         apply_to_structure(structure_A, new_P)
    #         writepdb("new_P_$(i).pdb", structure_A)
    #     else
    #         iterations_since_last_change += 1
    #     end

    #     # * Check if the error is less than the threshold.
    #     if error < ϵ
    #         apply_to_structure(structure_A, new_P)
    #         # * Write out the newly transformed PDB file.
    #         writepdb("new_P_$(i).pdb", structure_A)
    #         break
    #     end

    #     if iterations_since_last_change > 50
    #         # * If the error is not decreasing, perturb the system by introducing a random rotation.
    #         R = random_rotation()
    #         new_P = transpose(R * transpose(matrix_A) .+ t)
    #         matrix_A = new_P
    #         apply_to_structure(structure_A, new_P)
    #         writepdb("new_P_$(i).pdb", structure_A)
    #         iterations_since_last_change = 0
    #     end

    # end    


#     return 0

# end



# function kabsch_reject(x, y, max_iter::Int64, min_atoms::Int64)

#     # * truncate the matrices to the same size
#     x = x[1:min(size(x)[1], size(y)[1]), :]
#     y = y[1:min(size(x)[1], size(y)[1]), :]
    

#     n = size(x, 2)
#     best_rmsd = Inf
#     best_idx = collect(1:n)
#     for iter = 1:max_iter
#         # Randomly select a subset of atoms
#         idx = randperm(n)[1:min(n, convert(Int64, max(min_atoms, ceil(n / 2))))]

#         # Calculate the optimal rotation matrix and aligned coordinates
#         x_sub = x[:, idx]
#         y_sub = y[:, idx]
#         C = x_sub * y_sub'
#         U, S, V = svd(C)
#         d = det(U * V')
#         if d < 0
#             V[:, end] = -V[:, end]
#             d = -d
#         end
#         R = V * U'
#         y_aligned = R * y

#         # Calculate the RMSD for the subset of atoms
#         rmsd_sub = rmsd(x_sub, y_sub)

#         # If the RMSD for the subset is better than the best so far, update the best RMSD and index set
#         if rmsd_sub < best_rmsd
#             best_rmsd = rmsd_sub
#             best_idx = idx
#         end
#     end

#     # Align the structures using the best subset of atoms
#     x_sub = x[:, best_idx]
#     y_sub = y[:, best_idx]
#     C = x_sub * y_sub'
#     U, S, V = svd(C)
#     d = det(U * V')
#     if d < 0
#         V[:, end] = -V[:, end]
#         d = -d
#     end
#     R = V * U'
#     y_aligned = R * y

#     # Calculate the RMSD for the best subset of atoms
#     rmsd_best = rmsd(x_sub, y_sub)

#     return y_aligned, rmsd_best, best_idx
# end





# ****************************************************************************


# * Main * #

if abspath(PROGRAM_FILE) == @__FILE__

    CWD::String = pwd()
    file_A::String= String(ARGS[1])
    file_B::String = String(ARGS[2])
    A = read(file_A, PDB)
    B = read(file_B, PDB)
    pdb_coords_A = get_coords(file_A)
    pdb_coords_B = get_coords(file_B)
    iterations = parse(Int64, ARGS[3])

    kabsch_alignment(pdb_coords_A, pdb_coords_B, iterations, A, B)
    # kabsch_reject(pdb_coords_A, pdb_coords_B, iterations, 10)

end

