using BioStructures
using LinearAlgebra
using Statistics
using CoordinateTransformations, Rotations, StaticArrays






# ****************************************************************************
# * Point cloud structs -> https://docs.julialang.org/en/v1/manual/constructors/

struct PointCloud
    # points::Array{Tuple{Float64, Float64, Float64}, 1}
    points::Array{ Vector{Float64}, 1 } # * 2D array of 1D Vectors
    # *  Initializes a new instance of the struct PointCloud and sets the field "points" to an empty array.
    PointCloud() = new([])
end

function add_point(p::PointCloud, x::Float64, y::Float64, z::Float64)
    push!(p.points, [x, y, z])
end

function AddCoordsToPointCloud(pc::PointCloud, coords::Array{Float64, 2})
    for i in 1:size(coords, 1)
        add_point(pc, coords[i, 1], coords[i, 2], coords[i, 3])
    end
end

function print_pc(self::PointCloud)
    println("x\ty\tz")
    for coord_t in self.points
        println("$(coord_t[1])\t$(coord_t[2])\t$(coord_t[3])")
    end
end

# ****************************************************************************

function shortest_length(x::PointCloud, y::PointCloud)
    return min( length(x.points),length(y.points) )
end

# Function to extract coordinates from PDB file
function get_coords(pdbfile::String)
    # * Read in the structure and collect all of the atoms
    structure = read(pdbfile, PDB)
    atoms = collectatoms(structure)
    coordinates = Array{Float64}(undef, length(atoms), 3)

    # * Now, obtain the coordinates for each atom
    for (i, atom) in enumerate(atoms)
        coordinates[i, :] = atom.coords
    end
    
    return coordinates
end

# function rotation_trans_mat(rows::Int64, cols::Int64)
#     R = Matrix{Float64}(I, rows, cols)
#     t = Vector{Float64}(undef, rows)
#     return R, t
# end



function calculate_PC_centroid(pc::Matrix{Vector{Float64}})::Vector{Float64}

    # * Initialize the centroid to the zero vector
    centroid = zeros(Float64, 3)

    # * Calculate the centroid
    for coord_t in pc
        centroid += coord_t
    end
    centroid = centroid ./ length(pc)

    return centroid

end

function min_pts_AB_from_dist(matrix_in, min_dist)
    # min_dist_pts, true_min_pts  = [], []
    # for row in eachrow(matrix_in)
    #     # println("Last entry of the current row is $(row[3])")
    #     if row[2] == length(B.points)
    #         C = minimum(min_dist_pts)
    #         min_val = findall(isequal(minimum(min_dist_pts)), matrix_in)
    #         display(min_val)
    #         break
    #     else
    #         push!(min_dist_pts, row[3])
    #     end

    # end

    # display(min_dist_pts)

    # return min_dist_pts, true_min_pts


    for row in matrix_in
        if row[3] == min_dist
            point_A = row[1]
            point_B = row[2]
            break
        end
    end
   
    return point_A, point_B

end




function pt_correspondence(source::PointCloud, target::PointCloud)


    n = length(source.points)
    m = length(target.points)
    dist_mat = Matrix{Float64}(undef, n*m, 3)
    # variable_core_mat = Matrix{Float64}{undef, n, 3}

    i = 1
    for (idx_1, atom_A) in enumerate(source.points)
        for(idx_2, atom_B) in enumerate(target.points)
            x₁, x₂ = atom_A[1,1], atom_B[1,1]
            y₁, y₂ = atom_A[2,1], atom_B[2,1]
            z₁, z₂ = atom_A[3,1], atom_B[3,1]
            distance = sqrt( (x₁ - x₂)^2 + (y₁ - y₂)^2 + (z₁ - z₂)^2 )
            dist_mat[i, :] = [idx_1, idx_2, distance]
            i += 1
        end

        ##-- Filter the distance matrix find the closest point in A to B.
        ##-- Then, update the variable core matrix with the closest point.

    end


    # ##-- Need to find which point in A is closest to B.
    # y = findall( x -> x[1] == length(point_cloud_A.points), min_dist_pts )
    # first_y = first(y)[1]     #-- Isolating the first component of the first instance of the length in the cartesian object.
    # last_y = last(y)[1] 


    for i in 1:max(n,m)
        # println(i)
        y = findall( x -> x[1] == i, dist_mat )
        first_y = first(y)[1]
        last_y = last(y)[1]
        point_distances = dist_mat[first_y[1]:last_y[1], : ]
        display(point_distances)

        # min_third = minimum(row[3] for row in point_distances)
        # println(min_third)
    
        # println("$(first_y)\t$(last_y)")
        # println(dist_mat[i, :])
    end





    # min_dist_pts, true_min_pts  = [], []

    # for row in eachrow(dist_mat)
    #     # println("Last entry of the current row is $(row[3])")
    #     if row[2] == length(target.points)
    #         C = minimum(min_dist_pts)
    #         min_val = findall(isequal(minimum(min_dist_pts)), dist_mat)
    #         display(min_val)
    #         break
    #     else
    #         push!(min_dist_pts, row[3])
    #     end

    # end

    # display(min_dist_pts)

        # i += 1
        
        #     for (first_idx, vector) in enumerate(dist_mat)
        #         if vector[3] == minimum([ vector[3] for vector in dist_mat ])
        #             # println("$(vector[1])\t$(vector[2])\t$(vector[3])")
        #             push!(min_dist_mat, [vector[1], vector[2], vector[3]])
        #         end
        #     end
        
        # end
        

    
    
    # display(dist_mat)

    
    # return dist_mat
    # return 0 

    return dist_mat

end



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
    
    # * Adding coordinates to point each respective point cloud
    point_cloud_A = PointCloud()
    point_cloud_B = PointCloud()
    AddCoordsToPointCloud(point_cloud_A, pdb_coords_A)
    AddCoordsToPointCloud(point_cloud_B, pdb_coords_B)

    min_dist_pts = pt_correspondence(point_cloud_A, point_cloud_B)
    min_prot_len = min( length(point_cloud_A.points), length(point_cloud_B.points) )
    # display(min_dist_pts)

    ##-- Proof of concept: Isolating a set of point comparisons between
    ##-- the two point clouds, where the left hand side is the longest one,
    ##-- and the shortest one is in the middle column.
    # y = findall( x -> x[1] == length(point_cloud_A.points), min_dist_pts )
    # first_y = first(y)[1]     #-- Isolating the first component of the first instance of the length in the cartesian object.
    # last_y = last(y)[1]       #-- Isolating the first component of the legnth in the cartesian object.
    # display(min_dist_pts[first_y[1]:last_y[1], : ])








    
    # println("$(first_y)\t$(last_y)")



    # display(min_dist_pts[x])

    # x = min_dist_pts[findall( x -> x[1] == length(point_cloud_A.points), min_dist_pts), :]
    # display(x)

    # x = min_dist_pts[findall(x->x[1] ==1), :]
    # println(x)

    # A[findall(x -> x[1] == 1), :]

    # for i in 1:min_prot_len
        
    # end






    # display(min_dist_pts)
    # display(min_dist_pts[20:30, 1:3])


    


    
    # min_val, min_idx = findmin(distances[:,3])
    # println("$(min_idx)\t$(min_val)")
    
    # min_value = minimum([ vector[3] for vector in distances ])
    
    

    # println(A_coord_1)
    # for atom in point_cloud_B.points
    #     # distance = sqrt( sum( A_coord_1 - atom ).^2 )
    #     # println(distance)
    #     println(atom)
    # end














    # # * Initialize the transformation matrix to the identity matrix
    # R, t = rotation_trans_mat(3, 3)
    # display(R)
    # display(t)

  
    # # # * Calculate the centroid of the point clouds
    # pc_A_centroid = calculate_PC_centroid(point_cloud_A)
    # display(pc_A_centroid)


    # println(point_cloud_A)
    # print(pdb_coords_A[1, :])


    # # * Convert generator object to list
    # ref_structure_atoms = Vector(ref_structure.get_atoms())
    # sample_structure_atoms = Vector(sample_structure.get_atoms())

    # display(ref_structure_atoms)

    # pyplot()
    # p = scatter(x, y, z, markersize=0.5, markercolor="blue")
    # display(p)
    # gui()

end










# using LinearAlgebra
# using Random

# function icp(source::Array{Float64, 2}, target::Array{Float64, 2})
#     # Initialize the transformation matrix
#     T = eye(4)
#     # The termination criteria
#     epsilon = 1e-6
#     max_iter = 100
#     # Start the iteration
#     for i = 1:max_iter
#         # Find the closest points
#         distances = sum((source .- target').^2, dims=1)
#         closest_indices = argmin(distances, dims=1)
#         closest_points = target[:, closest_indices]
#         # Compute the centroids of the two point clouds
#         source_centroid = mean(source, dims=2)
#         closest_centroid = mean(closest_points, dims=2)
#         # Compute the difference between the centroids
#         translation = closest_centroid - source_centroid
#         # Compute the rotation matrix
#         source_centered = source .- source_centroid
#         closest_centered = closest_points .- closest_centroid
#         W = source_centered * closest_centered'
#         U, S, V = svd(W)
#         R = V * U'
#         # Update the transformation matrix
#         T[1:3, 1:3] = R
#         T[1:3, 4] = translation
#         # Transform the source point cloud
#         source = T * hcat(source, ones(size(source, 2)))
#         # Check the termination criteria
#         if norm(translation) < epsilon
#             break
#         end
#     end
#     return T
# end

# # Example usage
# Random.seed!(1)
# source = rand(3, 100)
# target = rand(3, 100)
# T = icp(source, target)






# function ICP(x::PointCloud, y::PointCloud)
#     # * Initialize the transformation matrix to the identity matrix
#     T = Matrix{Float64}(I, 3, 3)
#     # # * Initialize the translation vector to the zero vector
#     # t = zeros(3)
#     # # * Initialize the error vector to a vector of nans
#     # e = fill(NaN, 3)
#     # # * Initialize the error to a large number
#     # e_old = 1e6
#     # # * Initialize the error to a small number
#     # e_new = 1e-6
#     # # * Initialize the maximum number of iterations
#     # max_iterations = 100
#     # # * Initialize the iteration counter
#     # iteration = 0
#     # # * Initialize the convergence criteria
#     # convergence = 1e-6

#     # # * While the error is decreasing and the maximum number of iterations has not been reached
#     # while (e_new < e_old) && (iteration < max_iterations)
#     #     # * Update the old error
#     #     e_old = e_new
#     #     # * Increment the iteration counter
#     #     iteration += 1
#     #     # * Find the nearest neighbor in y for each point in x
#     #     # * Transform x
#     #     # * Update T and t
#     #     # * Compute the error
#     # end

# end


# # ! I think the error is occuring because I am trying to pass in a 2D array of 1D vectors,
# # ! but the function is expecting a 2D array of 3D vectors. I can fix this by changing the
# # ! struct PointCloud to have a 2D array of 3D vectors, but I am not sure if that is the
# # ! best way to do it. I will try it and see if it works. Or, I could just change the
# # ! function to accept a 2D array of 1D vectors, but I am not sure if that is the best course
# # ! of action either. I will try it and see if it works.
# function ICP(source::Matrix{Vector{Float64}}, target::Matrix{Vector{Float64}})
#     # * Initialize the transformation matrix to the identity matrix, unless otherwise specified
#     T = Matrix{Float64}(I, 4, 4)

#     # * Define Termination Criteria
#     ϵ::Float64 = 1e-6
#     max_iterations::Int64 = 50

#     println(size(source))
#     println(size(target))
    

#     # # * Start the iteration
#     # for i = 1:max_iterations
#     #     # Find the closest points
#     #     distances = sum((source .- target').^2, dims=1)
#     #     closest_indices = argmin(distances, dims=1)
#     #     closest_points = target[:, closest_indices]
#     # end

#     return T
# end

