abstract type  AbstractIOSetup end

mutable struct IOSetup <: AbstractIOSetup
    identifier :: Dict{Symbol, Vector{Int64}}     # Identifier
    indices :: Dict{Symbol, Vector{Int64}}   # Cell indices
    character :: Dict{Symbol, Vector{Real}}    # Configuration / characterization
end

function initIOSetup(geometry :: AbstractGeometricalObject)

    position_symbols = getboundarypositions( geometry )
    
    identifier = Dict{Symbol, Vector{Int64}}()
    indices = Dict{Symbol, Vector{Int64}}()
    character = Dict{Symbol, Vector{Real}}()

    for pos in position_symbols
        identifier[pos] = Int64[]
        indices[pos] = Int64[]
        character[pos] = Real[]
    end

    return IOSetup(identifier,indices,character)
end

function setIOSetup!(iosetup :: IOSetup, heatrod :: HeatRod, identifier_index :: Int64, config :: Real,  orientation :: Symbol)
    
    num_cells = heatrod.heatcells

    iosetup.identifier[orientation] = [identifier_index]
    iosetup.character[orientation] = [config]
    
    if orientation == :west
        iosetup.indices[orientation] = [1];
    elseif orientation == :east
        iosetup.indices[orientation]  = [num_cells];
    end

end


function setIOSetup!(iosetup :: IOSetup, heatplate :: HeatPlate, num_partitions :: Integer, config :: RadialConfiguration,  orientation :: Symbol; start_index = 1 :: Integer)
    
    partition_array = collect(start_index: 1: num_partitions+start_index-1)
    
    config_array = RadialConfiguration[]

    for i = 1 : num_partitions
        config_array = vcat(config_array, config)
    end

    setIOSetup!(iosetup , heatplate, partition_array, config_array,  orientation)

end

function setIOSetup!(iosetup :: IOSetup, heatplate :: HeatPlate, partition :: Vector{T} where T <: Integer, config :: Array{RadialConfiguration,1},  orientation :: Symbol)
    
    if size(partition) != size(config)
        error("Size of partition and configuration array have to be equal!")
    end


    if orientation == :south || orientation == :north
        step        = heatplate.sampling[1]
        distance    = heatplate.dimension[1]
        dim_index = 1
    elseif  orientation == :east || orientation == :west
        step     = heatplate.sampling[2]
        distance = heatplate.dimension[2]
        dim_index = 2
    else
        error("Orientation $(orientation) not available!")
    end
   
    num_partitions = length(partition)
    cellindices = getindices(heatplate, cellPosition = orientation) 
    centralpoints = findcenterpoints(heatplate, num_partitions, orientation) 

    if length(cellindices) < num_partitions 
        error("Number of partitions higher than number of cells! Partitions: $(num_partitions), Cells: $(length(cellindices)). \n You may change the size of partition and configuration array.")    end


    partitionrange = length(cellindices)/num_partitions;

   
    span = step/2 : step : distance - (step/2)
    stop_idx = 0;

    indices = Int64[]
    chars   = Real[]
    idents  = Int64[]

    for i = 1 : num_partitions

        if partition[i] > 0
            partitionindices = round(Int, (i-1) * partitionrange)+1 : 1 : round(Int, i * partitionrange);
            start_idx = partitionindices[1]
            stop_idx = partitionindices[end]
            
            xpos = centralpoints[1][i]
            ypos = centralpoints[2][i]
            zpos = centralpoints[3][i]
            point = (xpos, ypos, zpos);
            config[i].center = point
            
            spatial_character  = characterize(span[start_idx],span[stop_idx], step, config[i], dim=dim_index)

            identifier_index = partition[i]*ones(Int64, length(spatial_character))


            indices = vcat(indices, cellindices[partitionindices])
            chars = vcat(chars, spatial_character);
            idents = vcat(idents, identifier_index)
        end
    end

    iosetup.identifier[orientation] = idents;
    iosetup.character[orientation] = chars;
    iosetup.indices[orientation] = indices;
end



# HeatCuboid
function setIOSetup!(iosetup :: IOSetup, heatcuboid :: HeatCuboid, num_partitions ::  Tuple{T,T} , config :: RadialConfiguration,  orientation :: Symbol; start_index = 1 :: Integer) where T <: Integer
    
    num_actuators_row, num_actuators_col = num_partitions

    config_array    = Array{RadialConfiguration}(undef,num_actuators_row, num_actuators_col)
    partition_array = zeros(Integer,num_actuators_row, num_actuators_col)
    
    for iy = 1 : num_actuators_col
        for ix = 1 : num_actuators_row
            config_array[ix,iy] = config
            partition_array[ix,iy] = ix + (iy-1)*num_actuators_row  + (start_index-1)
        end
    end


    setIOSetup!(iosetup , heatcuboid, partition_array, config_array,  orientation)
end



# HeatCuboid
function setIOSetup!(iosetup :: IOSetup, heatcuboid :: HeatCuboid, partition :: Array{T,2} where T <: Integer, config :: Array{RadialConfiguration,2},  orientation :: Symbol)
    

    step = zeros(2)
    distance = zeros(2)
   
    if size(partition) != size(config)
        error("Size of partition and configuration array have to be equal!")
    end

    if orientation == :south || orientation == :north
        step[1] = heatcuboid.sampling[1]
        step[2] = heatcuboid.sampling[3]
        distance[1] = heatcuboid.dimension[1]
        distance[2] = heatcuboid.dimension[3]
        
        dim_index = (1,3)
    elseif  orientation == :east || orientation == :west
        step[1] = heatcuboid.sampling[2]
        step[2] = heatcuboid.sampling[3]

        distance[1] = heatcuboid.dimension[2]
        distance[2] = heatcuboid.dimension[3]
        dim_index = (2,3)

    elseif orientation == :underside || orientation == :topside
        step[1] = heatcuboid.sampling[1]
        step[2] = heatcuboid.sampling[2]
        distance[1] = heatcuboid.dimension[1]
        distance[2] = heatcuboid.dimension[2]
        
        dim_index = (1,2)
    else
        error("Orientation $(orient) not available!")
    end


    num_parts_rows, num_parts_cols = num_partitions = size(partition)
    centralpoints = findcenterpoints(heatcuboid, num_partitions, orientation) 
   

    num_boundarycells = getnumofboundarycells(heatcuboid, orientation)
    num_cells_row = num_boundarycells[1]
    num_cells_col = num_boundarycells[2]

    if num_cells_row < num_parts_rows 
        error("Number of row partitions higher than number of row cells! Partitions: $(num_parts_rows), Cells: $(num_cells_row). \n You may change the size of partition and configuration array.")
    elseif num_cells_col < num_parts_cols
        error("Number of column partitions higher than number of column cells! Partitions: $(num_parts_cols), Cells: $(num_cells_col). \n You may change the size of partition and configuration array.")
    end



    cellindices = getindices(heatcuboid, cellPosition = orientation) 
    cellindices = reshape(cellindices, num_cells_row, num_cells_col)



    partitionrange = zeros(2)
    partitionrange[1] = num_cells_row/num_parts_rows;
    partitionrange[2] = num_cells_col/num_parts_cols;

   
    span_row = step[1]/2 : step[1] : distance[1] - (step[1]/2)
    span_col = step[2]/2 : step[2] : distance[2] - (step[2]/2)
    
    
    stop_idx = 0;

    indices = Int64[]
    chars   = Real[]
    idents  = Int64[]


    for iy = 1 : num_parts_cols
        for ix = 1 : num_parts_rows
            if partition[ix,iy] > 0
                idx = ix + (iy-1)*num_parts_rows; # global index

                indices_row = round(Int, (ix-1) * partitionrange[1])+1 : 1 : round(Int, ix * partitionrange[1]);
                indices_col = round(Int, (iy-1) * partitionrange[2])+1 : 1 : round(Int, iy * partitionrange[2]);
            
                xpos = centralpoints[1][idx]
                ypos = centralpoints[2][idx]
                zpos = centralpoints[3][idx]
                point = (xpos, ypos, zpos);
                config[ix, iy].center = point

                span_start = (span_row[indices_row[1]], span_col[indices_col[1]])
                span_stop = (span_row[indices_row[end]], span_col[indices_col[end]])
                            
                spatial_character  = characterize(span_start, span_stop, (step[1], step[2]), config[ix,iy], dim=dim_index)
               
                num_cells    = length(spatial_character) # number of cells per partition
                identifier_index = partition[ix,iy]*ones(Int64, num_cells)

                indices = vcat(indices, reshape(cellindices[indices_row,indices_col],num_cells,1)[:,1])
                chars = vcat(chars, reshape(spatial_character, num_cells, 1)[:,1]);
                idents = vcat(idents, identifier_index)
            end
        end
    end

    iosetup.identifier[orientation] = idents;
    iosetup.character[orientation] = chars;
    iosetup.indices[orientation] = indices;

end

#=
    Returns a 2D spatial distribution of characterization and input indices
=#
function checkIOSetup2D(iosetup :: IOSetup, heatcuboid :: HeatCuboid, orientation :: Symbol)

    boundary_indices = getindices( heatcuboid ; cellPosition = orientation )
    num_boundarycells = getnumofboundarycells(heatcuboid, orientation)
    Nx = num_boundarycells[1]
    Ny = num_boundarycells[2]

    character = zeros(Nx,Ny);
    indices = zeros(Int64, Nx,Ny);
    
    for iy = 1 : Ny
        for ix = 1 : Nx
            i = ix + (iy-1)*Nx  
            cellindex = boundary_indices[i] 
            char, idx = getActuation(iosetup, cellindex, orientation)

            character[ix,iy] = char
            indices[ix,iy]   = idx
        end
    end

    return (character, indices)
end



function findcenterpoints(heatplate :: HeatPlate, num_elements :: Integer, position :: Symbol)
    
    partitionrange = nothing

    len  = heatplate.dimension[1]
    wid  = heatplate.dimension[2]

     
    if position == :south
        partitionrange = len / num_elements
        xpos = collect(partitionrange/2 : partitionrange : partitionrange*(2*num_elements - 1)/2)
        ypos = zeros(num_elements);

    elseif position == :north
        partitionrange = len / num_elements
        xpos = collect(partitionrange/2 : partitionrange : partitionrange*(2*num_elements - 1)/2)
        ypos = ones(num_elements) * wid;

    elseif position == :west
        partitionrange = wid/ num_elements
        xpos = zeros(num_elements)
        ypos = collect(partitionrange/2 : partitionrange : partitionrange*(2*num_elements - 1)/2)
    
    elseif position == :east
        partitionrange = wid/ num_elements
        xpos = ones(num_elements)*len
        ypos = collect(partitionrange/2 : partitionrange : partitionrange*(2*num_elements - 1)/2)
    
    else
        error("CellPosition $(position) is not available!")
    end
    

    zpos = zeros(num_elements);
    centralpoints = (xpos, ypos, zpos)

    return centralpoints
end





function findcenterpoints(heatcuboid :: HeatCuboid, num_elements :: Tuple{T,T}, position :: Symbol) where {T <: Integer}
  
    partitionrange = zeros(2);
    num_elements_total = num_elements[1]*num_elements[2];

    len  = heatcuboid.dimension[1] # length
    wid  = heatcuboid.dimension[2] # width
    hei  = heatcuboid.dimension[3] # heigh
     
    if position == :south
        partitionrange[1] = len / num_elements[1]
        partitionrange[2] = hei / num_elements[2]

        xpos = collect(partitionrange[1]/2 : partitionrange[1] : partitionrange[1]*(2*num_elements[1] - 1)/2)
        xpos = repeat(xpos, num_elements[2])

        ypos = zeros(num_elements_total);

        zpos = collect(partitionrange[2]/2 : partitionrange[2] : partitionrange[2]*(2*num_elements[2] - 1)/2)
        zpos = repeat(zpos,1,num_elements[1])
        zpos = reshape(zpos', length(zpos), 1)

    elseif position == :north
        partitionrange[1] = len / num_elements[1]
        partitionrange[2] = hei / num_elements[2]
        xpos = collect(partitionrange[1]/2 : partitionrange[1] : partitionrange[1]*(2*num_elements[1] - 1)/2)
        xpos = repeat(xpos, num_elements[2])

        ypos = ones(num_elements_total) * wid;
        
        zpos = collect(partitionrange[2]/2 : partitionrange[2] : partitionrange[2]*(2*num_elements[2] - 1)/2)
        zpos = repeat(zpos,1,num_elements[1])
        zpos = reshape(zpos', length(zpos), 1)

    elseif position == :west
        partitionrange[1] = wid/ num_elements[1]
        partitionrange[2] = hei / num_elements[2]

        xpos = zeros(num_elements_total)
        
        ypos = collect(partitionrange[1]/2 : partitionrange[1] : partitionrange[1]*(2*num_elements[1] - 1)/2)
        ypos = repeat(ypos, num_elements[2])

        zpos = collect(partitionrange[2]/2 : partitionrange[2] : partitionrange[2]*(2*num_elements[2] - 1)/2)
        zpos = repeat(zpos,1,num_elements[1])
        zpos = reshape(zpos', length(zpos), 1)

    elseif position == :east
        partitionrange[1] = wid/ num_elements[1]
        partitionrange[2] = hei / num_elements[2]

        xpos = ones(num_elements_total) * len

        ypos = collect(partitionrange[1]/2 : partitionrange[1] : partitionrange[1]*(2*num_elements[1] - 1)/2)
        ypos = repeat(ypos, num_elements[2])

        zpos = collect(partitionrange[2]/2 : partitionrange[2] : partitionrange[2]*(2*num_elements[2] - 1)/2)
        zpos = repeat(zpos,1,num_elements[1])
        zpos = reshape(zpos', length(zpos), 1)

        
    elseif  position == :underside
        partitionrange[1] = len / num_elements[1]
        partitionrange[2] = wid/ num_elements[2]

        xpos = collect(partitionrange[1]/2 : partitionrange[1] : partitionrange[1]*(2*num_elements[1] - 1)/2)
        xpos = repeat(xpos, num_elements[2])

        ypos = collect(partitionrange[2]/2 : partitionrange[2] : partitionrange[2]*(2*num_elements[2] - 1)/2)
        ypos = repeat(ypos,1,num_elements[1])
    
        ypos = reshape(ypos', length(ypos), 1)

        zpos = zeros(num_elements_total); 

    elseif position == :topside
        partitionrange[1] = len / num_elements[1]
        partitionrange[2] = wid/ num_elements[2]

        xpos = collect(partitionrange[1]/2 : partitionrange[1] : partitionrange[1]*(2*num_elements[1] - 1)/2)
        xpos = repeat(xpos, num_elements[2])

        ypos = collect(partitionrange[2]/2 : partitionrange[2] : partitionrange[2]*(2*num_elements[2] - 1)/2)
        ypos = repeat(ypos,1,num_elements[1])
        ypos = reshape(ypos', length(ypos), 1)

        zpos = ones(num_elements_total) * hei; 
        
    else
        error("CellPosition $(position) is not available!")
    end
    
    centralpoints = (xpos, ypos, zpos)

    return centralpoints
end




function getActuation(actuation :: IOSetup, cellindex :: Int64, orientation :: Symbol)

    index = findfirst(isequal(cellindex), actuation.indices[orientation])
    if isa(index, Nothing)
        b = 0;
        u_index = 1
        return (b, u_index)
    else
        b = actuation.character[orientation][index]
        u_index = actuation.identifier[orientation][index]
        return (b, u_index)
    end

end


function getSensing(sensor :: IOSetup, ioindex :: Int64, orientation :: Symbol)
   
    idents = sensor.identifier[orientation]

    ioindices = findall(isequal(ioindex), idents)
    cellindices = sensor.indices[orientation][ioindices]
    character   = sensor.character[orientation][ioindices]
    return (cellindices,character)
end


"""
Returns the weighted arithmetic mean of the measurement.
"""
function measureWAM(temperatures :: Vector{T1}, character :: Vector{T2}) where {T1 <: Real, T2 <: Real}
    return character' * temperatures / sum(character)
end

function measureWAM(temperatures :: Array{T1,2}, character :: Vector{T2}) where {T1 <: Real, T2 <: Real}
    return transpose(character' * temperatures / sum(character))
end