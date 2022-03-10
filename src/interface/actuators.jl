abstract type  AbstractActuation end
abstract type  AbstractBoundaryActuation <: AbstractActuation end



mutable struct HeatRodActuation <: AbstractBoundaryActuation
    indices_west :: Array{Integer,1}    # Contains 1 element
    indices_east :: Array{Integer,1}

    character_west :: Array{Real,1} # Contains 1 element
    character_east :: Array{Real,1}

    input_west :: Array{Integer,1} # Contains 1 element
    input_east :: Array{Integer,1}
end

mutable struct HeatPlateActuation <: AbstractBoundaryActuation
    indices_west :: Array{Integer,1}
    indices_east :: Array{Integer,1}
    indices_south :: Array{Integer,1}
    indices_north :: Array{Integer,1}

    character_west :: Array{Real,1}
    character_east :: Array{Real,1}
    character_south :: Array{Real,1}
    character_north :: Array{Real,1}

    input_west :: Array{Integer,1} # Contains the index of input signal
    input_east :: Array{Integer,1}
    input_south :: Array{Integer,1}
    input_north :: Array{Integer,1}
end


mutable struct HeatCuboidActuation <: AbstractBoundaryActuation
    indices_west :: Array{Integer,1}
    indices_east :: Array{Integer,1}
    indices_south :: Array{Integer,1}
    indices_north :: Array{Integer,1}
    indices_underside :: Array{Integer,1}
    indices_topside :: Array{Integer,1}

    character_west :: Array{Real,1}
    character_east :: Array{Real,1}
    character_south :: Array{Real,1}
    character_north :: Array{Real,1}
    character_underside :: Array{Real,1}
    character_topside :: Array{Real,1}

    input_west :: Array{Integer,1}
    input_east :: Array{Integer,1}
    input_south :: Array{Integer,1}
    input_north :: Array{Integer,1}
    input_underside :: Array{Integer,1}
    input_topside :: Array{Integer,1}
end



function initActuation(geometry :: AbstractGeometricalObject)
    

    if isa( geometry, HeatRod )
        iwest = Integer[]
        ieast = Integer[]

        charwest = Real[]
        chareast = Real[]

        inputwest = Integer[]
        inputeast = Integer[]

        indices = (iwest, ieast)
        characters = (charwest,chareast)
        inputs = (inputwest, inputeast)


        return HeatRodActuation(indices..., characters..., inputs...)

    elseif isa( geometry, HeatPlate )
        
        iwest = Integer[]
        ieast = Integer[] 
        isouth = Integer[]
        inorth = Integer[]
       
        charwest = Real[]
        chareast = Real[]
        charsouth = Real[]
        charnorth = Real[]
       
        inputwest = Integer[]
        inputeast = Integer[]
        inputsouth = Integer[]
        inputnorth = Integer[]
        
        indices = (iwest, ieast, isouth, inorth)
        characters = (charwest, chareast, charsouth, charnorth)
        inputs = (inputwest, inputeast, inputsouth, inputnorth)

        return HeatPlateActuation(indices..., characters..., inputs...)

    elseif isa( geometry, HeatCuboid )
          
        iwest = Integer[]
        ieast = Integer[]
        isouth = Integer[]
        inorth = Integer[]
        iunderside = Integer[]
        itopside = Integer[]

        charwest = Real[]
        chareast = Real[]
        charsouth = Real[]
        charnorth = Real[]
        charunderside = Real[]
        chartopside = Real[]
        
        inputwest = Integer[]
        inputeast = Integer[]
        inputsouth = Integer[]
        inputnorth = Integer[]
        inputunderside = Integer[]
        inputtopside = Integer[]

        indices = (iwest, ieast, isouth, inorth, iunderside, itopside)
        characters = (charwest, chareast, charsouth, charnorth, charunderside, chartopside)
        inputs = (inputwest, inputeast, inputsouth, inputnorth, inputunderside, inputtopside)

        return HeatCuboidActuation(indices..., characters..., inputs...)

    end

end

function setActuation!(actuation :: HeatRodActuation, heatrod :: HeatRod, input_index :: Integer, config :: Real,  orientation :: Symbol)
    
    num_cells = heatrod.heatcells

    if orientation == :west
        actuation.indices_west   = [1] 
        actuation.character_west = [config]   
        actuation.input_west     = [input_index]  
    elseif orientation == :east
        actuation.indices_east   = [num_cells] 
        actuation.character_east = [config]      
        actuation.input_east     = [input_index]  
    end

end


function setActuation!(actuation :: HeatPlateActuation, heatplate :: HeatPlate, num_partitions :: Integer, config :: RadialConfiguration,  orientation :: Symbol; start_index = 1 :: Integer)
    
    partition_array = collect(start_index: 1: num_partitions+start_index-1)
    
    config_array = RadialConfiguration[]

    for i = 1 : num_partitions
        config_array = vcat(config_array, config)
    end

    setActuation!(actuation , heatplate, partition_array, config_array,  orientation)

end

function setActuation!(actuation :: HeatPlateActuation, heatplate :: HeatPlate, partition :: Array{T,1} where T <: Integer, config :: Array{RadialConfiguration,1},  orientation :: Symbol)
    
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

    indices = Integer[]
    chars   = Real[]
    inputs  = Integer[]

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

            input_index = partition[i]*ones(Integer, length(spatial_character))


            indices = vcat(indices, cellindices[partitionindices])
            chars = vcat(chars, spatial_character);
            inputs = vcat(inputs, input_index)
        end
    end

    if orientation == :west
        actuation.indices_west   = indices 
        actuation.character_west = chars   
        actuation.input_west     = inputs  
    elseif orientation == :east
        actuation.indices_east   = indices 
        actuation.character_east = chars   
        actuation.input_east     = inputs
    elseif orientation == :south
        actuation.indices_south   = indices 
        actuation.character_south = chars   
        actuation.input_south     = inputs
    elseif orientation == :north
        actuation.indices_north   = indices 
        actuation.character_north = chars   
        actuation.input_north     = inputs
    end


end

function setActuation!(actuation :: HeatCuboidActuation, heatcuboid :: HeatCuboid, num_partitions ::  Tuple{T,T} , config :: RadialConfiguration,  orientation :: Symbol; start_index = 1 :: Integer) where T <: Integer
    
    num_actuators_row, num_actuators_col = num_partitions

    config_array    = Array{RadialConfiguration}(undef,num_actuators_row, num_actuators_col)
    partition_array = zeros(Integer,num_actuators_row, num_actuators_col)
    
    for iy = 1 : num_actuators_col
        for ix = 1 : num_actuators_row
            config_array[ix,iy] = config
            partition_array[ix,iy] = ix + (iy-1)*num_actuators_row  + (start_index-1)
        end
    end


    setActuation!(actuation , heatcuboid, partition_array, config_array,  orientation)
end



# HeatCuboid
function setActuation!(actuation :: HeatCuboidActuation, heatcuboid :: HeatCuboid, partition :: Array{T,2} where T <: Integer, config :: Array{RadialConfiguration,2},  orientation :: Symbol)
    

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

    indices = Integer[]
    chars   = Real[]
    inputs  = Integer[]


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
                input_index = partition[ix,iy]*ones(Integer, num_cells)

                indices = vcat(indices, reshape(cellindices[indices_row,indices_col],num_cells,1)[:,1])
                chars = vcat(chars, reshape(spatial_character, num_cells, 1)[:,1]);
                inputs = vcat(inputs, input_index)
            end
        end
    end


    if orientation == :west
        actuation.indices_west   = indices 
        actuation.character_west = chars   
        actuation.input_west     = inputs  
    elseif orientation == :east
        actuation.indices_east   = indices 
        actuation.character_east = chars   
        actuation.input_east     = inputs
    elseif orientation == :south
        actuation.indices_south   = indices 
        actuation.character_south = chars   
        actuation.input_south     = inputs
    elseif orientation == :north
        actuation.indices_north   = indices 
        actuation.character_north = chars   
        actuation.input_north     = inputs
    elseif orientation == :underside
        actuation.indices_underside   = indices 
        actuation.character_underside = chars   
        actuation.input_underside     = inputs
    elseif orientation == :topside
        actuation.indices_topside   = indices 
        actuation.character_topside = chars   
        actuation.input_topside     = inputs
    end


end




function getActuation(actuation :: AbstractBoundaryActuation, cellindex :: Integer, orientation :: Symbol)

    if orientation == :west
        index = findfirst(isequal(cellindex), actuation.indices_west)
        if isa(index, Nothing)
            b = 0;
            u_index = 1
            return (b, u_index)
        else
            b = actuation.character_west[index]
            u_index = actuation.input_west[index]
            return (b, u_index)
        end

    elseif orientation == :east
        index = findfirst(isequal(cellindex), actuation.indices_east)
        if isa(index, Nothing)
            b = 0;
            u_index = 1
            return (b, u_index)
        else
            b = actuation.character_east[index]
            u_index = actuation.input_east[index]
            return (b, u_index)
        end

    elseif orientation == :south
        index = findfirst(isequal(cellindex), actuation.indices_south)
        if isa(index, Nothing)
            b = 0;
            u_index = 1
            return (b, u_index)

        else
            b = actuation.character_south[index]
            u_index = actuation.input_south[index]
            return (b, u_index)
        end

    elseif orientation == :north
        index = findfirst(isequal(cellindex), actuation.indices_north)
        if isa(index, Nothing)
            b = 0;
            u_index = 1
            return (b, u_index)

        else
            b = actuation.character_north[index]
            u_index = actuation.input_north[index]
            return (b, u_index)
        end

    elseif orientation == :underside
        index = findfirst(isequal(cellindex), actuation.indices_underside)
        if isa(index, Nothing)
            b = 0;
            u_index = 1
            return (b, u_index)

        else
            b = actuation.character_underside[index]
            u_index = actuation.input_underside[index]
            return (b, u_index)
        end

    elseif orientation == :topside
        index = findfirst(isequal(cellindex), actuation.indices_topside)
        if isa(index, Nothing)
            b = 0;
            u_index = 1
            return (b, u_index)

        else
            b = actuation.character_topside[index]
            u_index = actuation.input_topside[index]
            return (b, u_index)
        end
    end

end






# Only for testing / not recommended yet
function sortActuation!(actuation :: HeatPlateActuation)
    
    if length(actuation.indices_west) > 0
        idx_perm = sortperm(actuation.indices_west)
        
        actuation.indices_west   = actuation.indices_west[idx_perm]
        actuation.character_west = actuation.character_west[idx_perm]   
        actuation.input_west     = actuation.input_west[idx_perm]  

    elseif length(actuation.indices_east) > 0
        idx_perm = sortperm(actuation.indices_east)

        actuation.indices_east   = actuation.indices_east[idx_perm] 
        actuation.character_east = actuation.character_east[idx_perm]   
        actuation.input_east     = actuation.input_east[idx_perm]

    elseif length(actuation.indices_south) > 0
        idx_perm = sortperm(actuation.indices_south)

        actuation.indices_south   = actuation.indices_south[idx_perm]
        actuation.character_south = actuation.character_south[idx_perm]  
        actuation.input_south     = actuation.input_south[idx_perm]

    elseif length(actuation.indices_north) > 0
        idx_perm = sortperm(actuation.indices_north)

        actuation.indices_north   = actuation.indices_north[idx_perm] 
        actuation.character_north = actuation.character_north[idx_perm]   
        actuation.input_north     = actuation.input_north[idx_perm]

    end

end

# Only for testing / not recommended yet
function sortActuation!(actuation :: HeatCuboidActuation)
    
    if length(actuation.indices_west) > 0
        idx_perm = sortperm(actuation.indices_west)
        
        actuation.indices_west   = actuation.indices_west[idx_perm]
        actuation.character_west = actuation.character_west[idx_perm]   
        actuation.input_west     = actuation.input_west[idx_perm]  

    elseif length(actuation.indices_east) > 0
        idx_perm = sortperm(actuation.indices_east)

        actuation.indices_east   = actuation.indices_east[idx_perm] 
        actuation.character_east = actuation.character_east[idx_perm]   
        actuation.input_east     = actuation.input_east[idx_perm]

    elseif length(actuation.indices_south) > 0
        idx_perm = sortperm(actuation.indices_south)

        actuation.indices_south   = actuation.indices_south[idx_perm]
        actuation.character_south = actuation.character_south[idx_perm]  
        actuation.input_south     = actuation.input_south[idx_perm]

    elseif length(actuation.indices_north) > 0
        idx_perm = sortperm(actuation.indices_north)

        actuation.indices_north   = actuation.indices_north[idx_perm] 
        actuation.character_north = actuation.character_north[idx_perm]   
        actuation.input_north     = actuation.input_north[idx_perm]

    elseif length(actuation.indices_underside) > 0
        idx_perm = sortperm(actuation.indices_underside)

        actuation.indices_underside   = actuation.indices_underside[idx_perm]  
        actuation.character_underside = actuation.character_underside[idx_perm]    
        actuation.input_underside     = actuation.input_underside[idx_perm] 

    elseif length(actuation.indices_topside) > 0
        idx_perm = sortperm(actuation.indices_topside)

        actuation.indices_topside   = actuation.indices_topside[idx_perm] 
        actuation.character_topside = actuation.character_topside[idx_perm]   
        actuation.input_topside     = actuation.input_topside[idx_perm]
    end

end


#=
    Returns a 2D spatial distribution of characterization and input indices
=#
function checkActuation2D(actuation :: HeatCuboidActuation, heatcuboid :: HeatCuboid, orientation :: Symbol)

    boundary_indices = getindices( heatcuboid ; cellPosition = orientation )
    num_boundarycells = getnumofboundarycells(heatcuboid, orientation)
    Nx = num_boundarycells[1]
    Ny = num_boundarycells[2]

    character = zeros(Nx,Ny);
    indices = zeros(Integer, Nx,Ny);
    
    for iy = 1 : Ny
        for ix = 1 : Nx
            i = ix + (iy-1)*Nx  
            cellindex = boundary_indices[i] 
            char, idx = getActuation(actuation, cellindex, orientation)

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
