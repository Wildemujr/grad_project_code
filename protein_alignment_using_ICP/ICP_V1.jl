# using GeometryTypes
# using Plots
# using PyCall
# @pyimport Bio.PDB as PDB
# @pyimport numpy as np
using BioStructures




# * Structs
# struct PointCloud
#     x::Array{Float64, 2}
#     y::Array{Float64, 2}
#     z::Array{Float64, 2}
# end

# mutable struct PointCloud
#     x::Float64
#     y::Float64
#     z::Float64
#     PointCloud(x,y,z) = new(x,y,z)
# end


# * Point cloud structs -> https://docs.julialang.org/en/v1/manual/constructors/
# * 
struct PointCloud
    points::Array{Tuple{Float64, Float64, Float64}, 1}
    # *  Initializes a new instance of the struct PointCloud 
    # *  and sets the field "points" to an empty array.
    PointCloud() = new([])
end

function add_point(p::PointCloud, x::Float64, y::Float64, z::Float64)
    push!(p.points, (x, y, z))
end

function print_pc(self::PointCloud)
    for coord_t in self.points
        println(coord_t)
    end
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








if abspath(PROGRAM_FILE) == @__FILE__
    
    # mutable struct PointCloud
    #     points::Array{Tuple{Float64, Float64, Float64}, 1}
    #     function PointCloud()
    #         points = Tuple{Float64, Float64, Float64}[]  # * Empty array of tuples
    #     end
    #     function add_point(p::PointCloud, x::Float64, y::Float64, z::Float64)
    #         push!(p.points, (x, y, z))
    #     end
    # end


    # struct PointCloud
    #     points::Array{Tuple{Float64,Float64,Float64}, 1}
    #     function PointCloud()
    #         points = Tuple{Float64,Float64,Float64}[]
    #     end
    #     function addpoint(p::PointCloud, x::Float64, y::Float64, z::Float64)
    #         push!(p.points, (x, y, z))
    #     end
    # end


    # struct PointCloud
    #     points::Array{Tuple{Float64,Float64,Float64}, 1}
    #     function PointCloud()
    #         points = Tuple{Float64,Float64,Float64}[]
    #     end
    #     function addpoint(p::PointCloud, x::Float64, y::Float64, z::Float64)
    #         push!(p.points, (x, y, z))
    #     end
    # end
    
    # struct PointCloud
    #     points::Array{Tuple{Float64,Float64,Float64}, 1}
    #     function PointCloud()
    #         points = Tuple{Float64,Float64,Float64}[]
    #         return points
    #     end
    #     function addpoint(p::PointCloud, x::Float64, y::Float64, z::Float64)
    #         push!(p.points, (x, y, z))
    #     end
    # end
    # pc1 = PointCloud()
    # pc1.addpoint(1.0, 2.0, 3.0)



    point_cloud_A = PointCloud()
    add_point(point_cloud_A, 1.0, 2.0, 3.0)
    add_point(point_cloud_A, 4.0, 5.0, 6.0)
    print_pc(point_cloud_A)

    # for i in pc.points
    #     println(i)
    # end


    CWD::String = pwd()
    file_A::String= String(ARGS[1])
    file_B::String = String(ARGS[2])
    A = read(file_A, PDB)
    B = read(file_B, PDB)

    pdb_coords_A = get_coords(file_A)
    pdb_coords_B = get_coords(file_B)



    # # * Use PDB Parser to read in the PDB files
    # parser = PDB.PDBParser()
    # ref_structure = parser.get_structure("reference", file_A)
    # sample_structure = parser.get_structure("sample", file_B)
    
    # # * Convert generator object to list
    # ref_structure_atoms = Vector(ref_structure.get_atoms())
    # sample_structure_atoms = Vector(sample_structure.get_atoms())

    # display(ref_structure_atoms)

    # pyplot()
    # p = scatter(x, y, z, markersize=0.5, markercolor="blue")
    # display(p)
    # gui()

end
