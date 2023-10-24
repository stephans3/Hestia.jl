abstract type AbstractGeometricalObject end
abstract type AbstractCubicObject <: AbstractGeometricalObject end


#############################################################

# 1D Heat Rod

#############################################################




"""
    HeatRod <: AbstractCubicObject
    
Model of an one dimensional rod for heat conduction

### Elements
`dimension` : length of the rod

`sampling` : spatial discretization 

`heatcells` : number of heatcells
"""

struct HeatRod <: AbstractCubicObject
    dimension :: Real         # Length of rod
    sampling  :: Real         # Spatial discretization: Δx
    heatcells :: Integer      # Number of heatcells 
end


"""
    HeatRod(rod_length :: Real,  heatcells :: Integer )

Returns a HeatRod model
"""
function HeatRod(rod_length :: Real,  heatcells :: Integer)
    if rod_length <= 0
        error("Length has to be greater than zero!");
        return nothing;
    end
    if heatcells < 1
        error("Number of elements of heatcells is less than one. Array of heat cells has to have at least one element.")
        return nothing;
    end

    dimension = rod_length;
    sampling = rod_length / heatcells

    return HeatRod(dimension, sampling, heatcells)
end




#############################################################

# 2D Heat plate

#############################################################

"""
    HeatPlate <: AbstractCubicObject
    
Model of a two dimensional plate for heat conduction

### Elements
`dimension` : tuple of length and width of the plate: (length, width)

`sampling` : tuple of spatial discretization: (Δx, Δy)

`heatcells` : Number of heatcells per direction: {Nx, Ny}
"""
struct HeatPlate <: AbstractCubicObject
    dimension :: Tuple{T1, T1}  where T1 <: Real    # {Length, Width} of plate
    sampling  :: Tuple{T2, T2}  where T2 <: Real    # Spatial discretization: {Δx, Δy}
    heatcells :: Tuple{T3, T3}  where T3 <: Integer # Number of heatcells per direction: {Nx, Ny}
end



"""
    HeatPlate(plate_length :: Real, plate_width ::  Real, Nx :: Integer, Ny :: Integer, heatcells :: Array{S,1} where S <: Real)

Returns a HeatPlate model
"""
function HeatPlate( plate_length :: Real, 
                    plate_width  ::  Real, 
                    Nx :: Integer, 
                    Ny :: Integer)

    if plate_length <= 0 || plate_width <= 0
        error("Length and width have to be greater than zero!");
        return nothing;
    end

    if Nx < 2 || Ny < 2
        error("Number of column and row elements have to be greater than two!");
        return nothing;
    end
 
    Δx = plate_length/Nx;
    Δy = plate_width/Ny;

    dimension = (plate_length, plate_width)     # {Length, Width} of plate
    sampling  = (Δx, Δy)                        # Spatial discretization: {Δx, Δy}
    heatcells = (Nx, Ny)
    
    return HeatPlate(dimension, sampling, heatcells)
end


#############################################################

# 3D Heat cuboid

#############################################################
"""
    HeatCuboid <: AbstractCubicObject
    
Model of a two dimensional plate for heat conduction

### Elements
`dimension` : tuple of length, width and heigth of the plate: (length, width, height)

`sampling` : tuple of spatial discretization: (Δx, Δy, Δz)

`heatcells` : number of heatcells in total
"""
struct HeatCuboid <: AbstractCubicObject
    dimension :: Tuple{T1, T1, T1}  where T1 <: Real    # {Length, Width, Height} of cuboid
    sampling  :: Tuple{T2, T2, T2}  where T2 <: Real    # Spatial discretization: {Δx, Δy, Δz}
    heatcells :: Tuple{T3, T3, T3}  where T3 <: Integer # Number of heatcells per direction: {Nx, Ny, Nz} 
end


function HeatCuboid(cuboid_length :: Real, 
                    cuboid_width  ::  Real, 
                    cuboid_height ::  Real,
                    Nx :: Integer, 
                    Ny :: Integer, 
                    Nz :: Integer)

    if cuboid_length <= 0 || cuboid_width <= 0 || cuboid_height <= 0
        error("Length and width have to be greater than zero!");
        return nothing;
    end

    if Nx < 2 || Ny < 2 || Nz < 2
        error("Number of column and row elements have to be greater than zero!");
        return nothing;
    end
   
    Δx = cuboid_length/Nx;
    Δy = cuboid_width/Ny;
    Δz = cuboid_height/Nz;

    dimension = (cuboid_length, cuboid_width, cuboid_height)    # {Length, Width} of plate
    sampling  = (Δx, Δy, Δz)                                    # Spatial discretization: {Δx, Δy}
    heatcells = (Nx, Ny, Nz)

    return HeatCuboid(dimension, sampling, heatcells)
end


function isposition(pos :: Symbol)
    pos_list = (:complete, :center, :lateral, :west, :east, :south, :north, :underside, :topside)

    if pos in pos_list
        return true
    else
        return false
    end
end

function getpositions( geometry :: AbstractCubicObject )
    if isa( geometry, HeatRod )
        return (:complete, :center, :lateral, :west, :east)

    elseif isa( geometry, HeatPlate )
        return (:complete, :center, :lateral, :west, :east, :south, :north)

    elseif isa( geometry, HeatCuboid )
        return (:complete, :center, :lateral, :west, :east, :south, :north, :underside, :topside)
    else
        return ()
    end
end

function getboundarypositions( geometry :: AbstractCubicObject )
    if isa( geometry, HeatRod )
        return (:west, :east)

    elseif isa( geometry, HeatPlate )
        return (:west, :east, :south, :north)

    elseif isa( geometry, HeatCuboid )
        return (:west, :east, :south, :north, :underside, :topside)
    else
        return ()
    end
end


"""
    getindices( heatrod :: HeatRod ; cellPosition::Symbol = :complete )

Returns all indices for a certain `cellPosition`

## Valid CellPositions

`:complete`

`:center`

`:west`

`:east`

"""
function getindices( heatrod :: HeatRod ; cellPosition::Symbol = :complete )

    listpositions = getpositions(heatrod)

    if isposition(cellPosition) == false || (cellPosition in listpositions) == false
        error("CellPosition=$cellPosition not valid position!")
    end

    Nx = heatrod.heatcells

    indices = Integer[]

    if cellPosition == :complete
        indices = collect(1 : heatrod.heatcells)
    end
    
    if cellPosition == :center
        if Nx < 3
            error("Indices can not be found. No center available.")
        end

        for i = 2 : Nx-1
            idx = i;
            push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :west)
        push!(indices, 1)
    end

    if cellPosition in (:lateral, :east)
        push!(indices, Nx)
    end

    return indices
end

"""
    getindices( heatPlate :: HeatPlate ; cellPosition::Symbol = :complete )

Returns all indices for a certain `cellPosition`

## Valid CellPositions

`:complete`

`:center`

`:west`

`:east`

`:south`

`:north`

"""
function getindices( heatPlate :: HeatPlate ; cellPosition::Symbol = :complete )

    listpositions = getpositions(heatPlate)

    if isposition(cellPosition) == false || (cellPosition in listpositions) == false
        error("CellPosition=$cellPosition not valid position!")
    end

    Nx, Ny  = heatPlate.heatcells

    indices = Integer[]


    if cellPosition == :complete
        indices = collect(1 : Nx*Ny)
    end
    
    if cellPosition == :center
        if Nx < 3 || Ny < 3
            error("Indices can not be found. No center available.")
        end

        for j = 2 : Ny-1 , i = 2 : Nx-1
            idx = i + Nx*(j-1);
            push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :west)
        for j = 1 : Ny
           idx = 1 + Nx*(j-1)
           push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :east)
        for j = 1 : Ny
            idx = Nx*j
            push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :south)
        idx = collect(1:Nx)
        indices = vcat(indices,idx)
    end

    if  cellPosition in (:lateral, :north)
        startindex = Nx*(Ny - 1) + 1
        finalindex = Nx*Ny
        idx = collect(startindex : 1 :  finalindex)
        indices = vcat(indices,idx)
    end


    return indices
end







function getindices( heatcuboid :: HeatCuboid ; cellPosition::Symbol = :complete )

    listpositions = getpositions(heatcuboid)

    if isposition(cellPosition) == false || (cellPosition in listpositions) == false
        error("CellPosition=$cellPosition not valid position!")
    end
    
    Nx, Ny, Nz = heatcuboid.heatcells

    indices = Integer[]


    if cellPosition == :complete
        indices = collect(1 : Nx*Ny*Nz)
    end
    
    if cellPosition == :center
        if Nx < 3 || Ny < 3 || Nz < 3
            error("Indices can not be found. No center available.")
        end

        for k = 2 : Nz-1, j = 2 : Ny-1 , i = 2 : Nx-1
            idx = i + Nx*(j-1) + Nx*Ny*(k-1);
            push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :west)
        for k = 1 : Nz, j = 1 : Ny
           idx = 1 + Nx*(j-1) + Nx*Ny*(k-1);
           push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :east)
        for k = 1 : Nz, j = 1 : Ny
            idx = Nx*j  + Nx*Ny*(k-1)
            push!(indices, idx)
        end
    end

    if cellPosition in (:lateral, :south)
        for k = 1 : Nz
            start_idx = 1 + (k-1)*Nx*Ny
            final_idx = start_idx + Nx -1
            idx = collect(start_idx : final_idx)
            indices = vcat(indices,idx)
        end
        
    end

    if  cellPosition in (:lateral, :north)
        for k = 1 : Nz
            start_idx = Nx*(Ny - 1) + 1 + (k-1)*Nx*Ny
            final_idx = start_idx + Nx -1
            idx = collect(start_idx : final_idx)
            indices = vcat(indices,idx)
        end
    end

    if  cellPosition in (:lateral, :underside)
        start_idx = 1
        final_idx = Nx * Ny
        idx = collect(start_idx : final_idx)
        indices = vcat(indices,idx)
    end

    if  cellPosition in (:lateral, :topside)
        start_idx = 1 + (Nz-1)*Nx*Ny
        final_idx = Nx * Ny * Nz
        idx = collect(start_idx : final_idx)
        indices = vcat(indices,idx)
    end


    return indices
end


function getnumofboundarycells(heatcuboid :: HeatCuboid, position :: Symbol)
    
    Nx, Ny, Nz = heatcuboid.heatcells

    if position in (:west, :east)
        return (Ny, Nz)
    
    elseif position in (:south, :north)
        return (Nx, Nz)
    
    elseif position in (:underside, :topside)
        return (Nx, Ny)
    
    else
        return (0,0)
    end
end


function getindexofposition( heatcuboid :: HeatCuboid, position :: Tuple{T,T,T}) where {T<:Integer}
    Nx, Ny, Nz = heatcuboid.heatcells

    j = position[1]
    m = position[2]
    k = position[3]

    idx = j + (m-1)*Nx + (k-1)*Nx*Ny

    return idx
end

function getpositionofindex(heatcuboid :: HeatCuboid, index :: Integer)

    Nx, Ny, Nz = heatcuboid.heatcells
    
    if index > Nx*Ny*Nz
        error("Index higher than total number of heatcells")
    end
  

    x_pos = (index-1) % Nx + 1
    z_pos = floor(Int, (index-1)/(Nx * Ny)) + 1
    y_pos = round(Int, (index - x_pos - (z_pos-1)*Ny*Nx)/Nx) + 1


    return (x_pos, y_pos, z_pos)

end
