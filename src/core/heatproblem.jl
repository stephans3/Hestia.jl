abstract type AbstractHeatProblem end

struct CubicHeatProblem <:AbstractHeatProblem
    geometry    :: AbstractCubicObject
    boundary    :: CubicBoundary
end

function diffusion!(dθ :: AbstractArray{<: Real}, 
    θ :: AbstractArray{<: Real},  
    heatrod :: HeatRod, 
    boundary :: CubicBoundary)

    
end


function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatrod :: HeatRod, 
                    property :: StaticIsotropic, 
                    boundary :: CubicBoundary) where {T1 <: Real, T2 <: Real}
    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density
    
    Nx = heatrod.heatcells
    Δx = heatrod.sampling
    
    diffusion_static_x!(dθ, θ, Nx, 1, 1, Δx, λ, c, ρ, boundary)
end



function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatplate :: HeatPlate, 
                    property :: StaticIsotropic, 
                    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}
    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny) = heatplate.heatcells 
    (Δx,Δy) = heatplate.sampling
        
    diffusion_static_x!(dθ, θ, Nx, Ny, 1 , Δx, λ, c, ρ, boundary)
    diffusion_static_y!(dθ, θ, Nx, Ny, 1, Δy, λ, c, ρ, boundary)
end

function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatcuboid :: HeatCuboid, 
                    property :: StaticIsotropic, 
                    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}
    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling
        
    diffusion_static_x!(dθ, θ, Nx, Ny, Nz, Δx, λ, c, ρ, boundary)
    diffusion_static_y!(dθ, θ, Nx, Ny, Nz, Δy, λ, c, ρ, boundary)
    diffusion_static_z!(dθ, θ, Nx, Ny, Nz, Δz, λ, c, ρ, boundary)
end



# Anisotropic heat conduction
function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatplate :: HeatPlate, 
    property :: StaticAnisotropic, 
    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}

    λx = property.λx # Thermal conductivity in x-direction
    λy = property.λy # ... y-direction

    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny) = heatplate.heatcells 
    (Δx,Δy) = heatplate.sampling

    diffusion_static_x!(dθ, θ, Nx, Ny, 1 , Δx, λx, c, ρ, boundary)
    diffusion_static_y!(dθ, θ, Nx, Ny, 1, Δy, λy, c, ρ, boundary)
end

function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatcuboid :: HeatCuboid, 
    property :: StaticAnisotropic, 
    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}

    λx = property.λx # Thermal conductivity in x-direction
    λy = property.λy # ... y-direction
    λz = property.λz # ... z-direction
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling

    diffusion_static_x!(dθ, θ, Nx, Ny, Nz, Δx, λx, c, ρ, boundary)
    diffusion_static_y!(dθ, θ, Nx, Ny, Nz, Δy, λy, c, ρ, boundary)
    diffusion_static_z!(dθ, θ, Nx, Ny, Nz, Δz, λz, c, ρ, boundary)
end



function diffusion_static_x!(dθ,θ,Nx, Ny, Nz, Δx, λ :: Real, c :: Real, ρ :: Real, boundary :: CubicBoundary) # in-place 
    α = λ/(c*ρ) # Diffusivity

    @inbounds for iz in 1:Nz      
        @inbounds for iy in 1 : Ny
            @inbounds for ix in 2 : Nx-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                
                dθ[i] = α * (θ[i-1] + θ[i+1] - 2*θ[i])/Δx^2 
            end
            i1 = (iz-1)*Nx*Ny + (iy-1)*Nx + 1      # West
            i2 = (iz-1)*Nx*Ny + (iy-1)*Nx + Nx     # East
            
            emission_west = getEmission(boundary, i1, :west)
            emission_east = getEmission(boundary, i2, :east)
    
            dθ[i1] = α * (θ[i1+1] - θ[i1]) / Δx^2 + emit(θ[i1], emission_west)/(c*ρ*Δx)
            dθ[i2] = α * (θ[i2-1] - θ[i2]) / Δx^2 + emit(θ[i2], emission_east)/(c*ρ*Δx)
        end
    end
    nothing 
end


function diffusion_static_y!(dθ,θ,Nx, Ny, Nz, Δy, λ :: Real, c :: Real, ρ :: Real, boundary :: CubicBoundary) # in-place
    α = λ/(c*ρ)
    
    @inbounds for iz in 1:Nz
        @inbounds for ix in 1 : Nx
            @inbounds for  iy in 2 : Ny-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                
                dθ[i] = dθ[i] + α * (θ[i-Nx] + θ[i+Nx] - 2*θ[i])/Δy^2 
    
            end
            i1 = (iz-1)*Nx*Ny + ix                 # South
            i2 = (iz-1)*Nx*Ny + Nx*(Ny-1) + ix     # North
            
            emission_south = getEmission(boundary, i1, :south)
            emission_north = getEmission(boundary, i2, :north)
    
            dθ[i1] = dθ[i1] + α * (θ[i1+Nx] - θ[i1]) / Δy^2 + emit(θ[i1], emission_south)/(c*ρ*Δy)
            dθ[i2] = dθ[i2] + α * (θ[i2-Nx] - θ[i2]) / Δy^2 + emit(θ[i2], emission_north)/(c*ρ*Δy)
        end
    end

    nothing 
end

function diffusion_static_z!(dθ,θ,Nx, Ny, Nz, Δz, λ :: Real, c :: Real, ρ :: Real, boundary :: CubicBoundary) # in-place
    α = λ/(c*ρ)
    
    @inbounds for ix in 1 : Nx
        @inbounds for  iy in 1 : Ny
            @inbounds for iz in 2:Nz-1

                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                
                dθ[i] = dθ[i] + α * (θ[i-Nx*Ny] + θ[i+Nx*Ny] - 2*θ[i])/Δz^2 
    
            end
            
            i1 = (iy-1)*Nx + ix                 # Underside
            i2 = (Nz-1)*Nx*Ny + (iy-1)*Nx + ix  # Topside
            
            emission_underside = getEmission(boundary, i1, :underside)
            emission_topside   = getEmission(boundary, i2, :topside)
    
            dθ[i1] = dθ[i1] + α * (θ[i1+Nx*Ny] - θ[i1]) / Δz^2 + emit(θ[i1], emission_underside)/(c*ρ*Δz)
            dθ[i2] = dθ[i2] + α * (θ[i2-Nx*Ny] - θ[i2]) / Δz^2 + emit(θ[i2], emission_topside  )/(c*ρ*Δz)
        end
    end

    nothing 
end





function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatrod :: HeatRod, 
                    property :: DynamicIsotropic, 
                    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}
    
    λ(x) = specifyproperty(x, property.λ) # thermal conductivity
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    Nx = heatrod.heatcells
    Δx = heatrod.sampling
    
    diffusion_dynamic_x!(dθ, θ, Nx, 1, 1, Δx, λ, c, ρ, boundary)
end


function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatplate :: HeatPlate, 
                    property :: DynamicIsotropic, 
                    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}
   
    λ(x) = specifyproperty(x, property.λ) # thermal conductivity
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny) = heatplate.heatcells 
    (Δx,Δy) = heatplate.sampling
        
    diffusion_dynamic_x!(dθ, θ, Nx, Ny, 1, Δx, λ, c, ρ, boundary)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, 1, Δy, λ, c, ρ, boundary)
end

function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatcuboid :: HeatCuboid, 
                    property :: DynamicIsotropic, 
                    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}
   
    λ(x) = specifyproperty(x, property.λ) # thermal conductivity
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling
        
    diffusion_dynamic_x!(dθ, θ, Nx, Ny, Nz, Δx, λ, c, ρ, boundary)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, Nz, Δy, λ, c, ρ, boundary)
    diffusion_dynamic_z!(dθ, θ, Nx, Ny, Nz, Δz, λ, c, ρ, boundary)
end




function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatplate :: HeatPlate, 
    property :: DynamicAnisotropic, 
    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}

    λx(x) = specifyproperty(x, property.λx) # thermal conductivity in x-direction
    λy(x) = specifyproperty(x, property.λy) # ... in y-direction

    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny) = heatplate.heatcells 
    (Δx,Δy) = heatplate.sampling

    diffusion_dynamic_x!(dθ, θ, Nx, Ny, 1, Δx, λx, c, ρ, boundary)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, 1, Δy, λy, c, ρ, boundary)
end


function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatcuboid :: HeatCuboid, 
    property :: DynamicAnisotropic, 
    boundary :: CubicBoundary ) where {T1 <: Real, T2 <: Real}

    λx(x) = specifyproperty(x, property.λx) # thermal conductivity in x-direction
    λy(x) = specifyproperty(x, property.λy) # ... in y-direction
    λz(x) = specifyproperty(x, property.λz) # ... in z-direction
    
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling

    diffusion_dynamic_x!(dθ, θ, Nx, Ny, Nz, Δx, λx, c, ρ, boundary)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, Nz, Δy, λy, c, ρ, boundary)
    diffusion_dynamic_z!(dθ, θ, Nx, Ny, Nz, Δz, λz, c, ρ, boundary)
end


function diffusion_dynamic_x!(dθ,θ,Nx, Ny, Nz, Δx, λ :: Function, c :: Function, ρ :: Function, boundary :: CubicBoundary ) # in-place
    
    @inbounds for iz in 1:Nz
        @inbounds for iy in 1 : Ny
            @inbounds for ix in 2 : Nx-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                x̅1 =  (θ[i-1] + θ[i])/2
                x̅2 =  (θ[i+1] + θ[i])/2
                
                dθ[i] = (λ(x̅1) * θ[i-1] + λ(x̅2) * θ[i+1] - (λ(x̅1) + λ(x̅2))*θ[i])/(Δx^2 * ρ(θ[i]) * c(θ[i]))
            end
            i1 = (iz-1)*Nx*Ny + (iy-1)*Nx + 1      # West
            i2 = (iz-1)*Nx*Ny + (iy-1)*Nx + Nx     # East
         
            x̅1 =  (θ[i1] + θ[i1+1])/2
            x̅2 =  (θ[i2] + θ[i2-1])/2
            
            emission_west = getEmission(boundary, i1, :west)
            emission_east = getEmission(boundary, i2, :east)
    
            dθ[i1] = λ(x̅1) * (θ[i1+1] - θ[i1]) / (Δx^2 * ρ(θ[i1]) * c(θ[i1])) + emit(θ[i1], emission_west)/(Δx * ρ(θ[i1]) * c(θ[i1]))
            dθ[i2] = λ(x̅2) * (θ[i2-1] - θ[i2]) / (Δx^2 * ρ(θ[i2]) * c(θ[i2])) + emit(θ[i2], emission_east)/(Δx * ρ(θ[i2]) * c(θ[i2]))
        end
    end

    

    nothing 
end



function diffusion_dynamic_y!(dθ,θ,Nx, Ny, Nz, Δy, λ :: Function, c :: Function, ρ :: Function, boundary :: CubicBoundary ) # in-place
    @inbounds for iz in 1:Nz
        @inbounds for ix in 1 : Nx
            @inbounds for  iy in 2 : Ny-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                x̅1 =  (θ[i-Nx] + θ[i])/2
                x̅2 =  (θ[i+Nx] + θ[i])/2
                dθ[i] = dθ[i] +  (λ(x̅1) * θ[i-Nx] + λ(x̅2) * θ[i+Nx] - (λ(x̅1) + λ(x̅2))*θ[i])/(Δy^2 * ρ(θ[i]) * c(θ[i]))
    
            end
            i1 = (iz-1)*Nx*Ny + ix                 # South
            i2 = (iz-1)*Nx*Ny + Nx*(Ny-1) + ix     # North
    
            x̅1 =  (θ[i1] + θ[i1+Nx])/2
            x̅2 =  (θ[i2] + θ[i2-Nx])/2
    
            emission_south = getEmission(boundary, i1, :south)
            emission_north = getEmission(boundary, i2, :north)
    
            dθ[i1] = dθ[i1] + λ(x̅1) * (θ[i1+Nx] - θ[i1]) / (Δy^2 * ρ(θ[i1]) * c(θ[i1])) + emit(θ[i1],emission_south)/ (Δy * ρ(θ[i1]) * c(θ[i1]))
            dθ[i2] = dθ[i2] + λ(x̅2) * (θ[i2-Nx] - θ[i2]) / (Δy^2 * ρ(θ[i2]) * c(θ[i2])) + emit(θ[i2],emission_north)/ (Δy * ρ(θ[i2]) * c(θ[i2]))
        end
    end   

    nothing 
end




function diffusion_dynamic_z!(dθ,θ,Nx, Ny, Nz, Δz,  λ :: Function, c :: Function, ρ :: Function, boundary :: CubicBoundary) # in-place
    
    @inbounds for ix in 1 : Nx
        @inbounds for  iy in 1 : Ny
            @inbounds for iz in 2:Nz-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                x̅1 =  (θ[i-Nx*Ny] + θ[i])/2
                x̅2 =  (θ[i+Nx*Ny] + θ[i])/2
                
                dθ[i] = dθ[i] +  (λ(x̅1) * θ[i-Nx*Ny] + λ(x̅2) * θ[i+Nx*Ny] - (λ(x̅1) + λ(x̅2))*θ[i])/(Δz^2 * ρ(θ[i]) * c(θ[i]))
            end
            
            i1 = (iy-1)*Nx + ix                 # Underside
            i2 = (Nz-1)*Nx*Ny + (iy-1)*Nx + ix  # Topside
            
            x̅1 =  (θ[i1] + θ[i1+Nx])/2
            x̅2 =  (θ[i2] + θ[i2-Nx])/2

            emission_underside = getEmission(boundary, i1, :underside)
            emission_topside   = getEmission(boundary, i2, :topside  )
    
            dθ[i1] = dθ[i1] + λ(x̅1) * (θ[i1+Nx*Ny] - θ[i1]) / (Δz^2 * ρ(θ[i1]) * c(θ[i1])) + emit(θ[i1],emission_underside)/ (Δz * ρ(θ[i1]) * c(θ[i1]))
            dθ[i2] = dθ[i2] + λ(x̅2) * (θ[i2-Nx*Ny] - θ[i2]) / (Δz^2 * ρ(θ[i2]) * c(θ[i2])) + emit(θ[i2],emission_topside  )/ (Δz * ρ(θ[i2]) * c(θ[i2]))
        end
    end

    nothing 
end



function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatrod :: HeatRod, 
                    property :: StaticIsotropic, 
                    boundary :: CubicBoundary,  
                    actuation :: AbstractIOSetup, 
                    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}

    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density
    
    Nx = heatrod.heatcells
    Δx = heatrod.sampling

    diffusion_static_x!(dθ, θ, Nx, 1, 1 , Δx, λ, c, ρ, boundary, actuation, input_signals)
end



function diffusion!(dθ :: AbstractArray{T1}, 
                        θ :: AbstractArray{T2},  
                        heatplate :: HeatPlate, 
                        property :: StaticIsotropic, 
                        boundary :: CubicBoundary,  
                        actuation :: AbstractIOSetup, 
                        input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}

    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny) = heatplate.heatcells 
    (Δx,Δy) = heatplate.sampling
        
    diffusion_static_x!(dθ, θ, Nx, Ny, 1 , Δx, λ, c, ρ, boundary, actuation, input_signals)
    diffusion_static_y!(dθ, θ, Nx, Ny, 1, Δy, λ, c, ρ, boundary, actuation, input_signals)
end

function diffusion!(dθ :: AbstractArray{T1}, 
                        θ :: AbstractArray{T2},  
                        heatcuboid :: HeatCuboid, 
                        property :: StaticIsotropic, 
                        boundary :: CubicBoundary,
                        actuation :: AbstractIOSetup, 
                        input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}
    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling
        
    diffusion_static_x!(dθ, θ, Nx, Ny, Nz, Δx, λ, c, ρ, boundary, actuation, input_signals)
    diffusion_static_y!(dθ, θ, Nx, Ny, Nz, Δy, λ, c, ρ, boundary, actuation, input_signals)
    diffusion_static_z!(dθ, θ, Nx, Ny, Nz, Δz, λ, c, ρ, boundary, actuation, input_signals)
end


function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatplate :: HeatPlate, 
    property :: StaticAnisotropic, 
    boundary :: CubicBoundary,  
    actuation :: AbstractIOSetup, 
    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}

    λx = property.λx # Thermal conductivity x-direction
    λy = property.λy # ... in y-direction
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny) = heatplate.heatcells 
    (Δx,Δy) = heatplate.sampling

    diffusion_static_x!(dθ, θ, Nx, Ny, 1 , Δx, λx, c, ρ, boundary, actuation, input_signals)
    diffusion_static_y!(dθ, θ, Nx, Ny, 1, Δy, λy, c, ρ, boundary, actuation, input_signals)
end

function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatcuboid :: HeatCuboid, 
    property :: StaticAnisotropic, 
    boundary :: CubicBoundary,
    actuation :: AbstractIOSetup, 
    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}

    λx = property.λx # Thermal conductivity x-direction
    λy = property.λy # ... in y-direction
    λz = property.λz # ... in z-direction
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling

    diffusion_static_x!(dθ, θ, Nx, Ny, Nz, Δx, λx, c, ρ, boundary, actuation, input_signals)
    diffusion_static_y!(dθ, θ, Nx, Ny, Nz, Δy, λy, c, ρ, boundary, actuation, input_signals)
    diffusion_static_z!(dθ, θ, Nx, Ny, Nz, Δz, λz, c, ρ, boundary, actuation, input_signals)
end




function diffusion_static_x!(dθ, θ, Nx, Ny, Nz, Δx, λ :: Real, c :: Real, ρ :: Real, boundary :: CubicBoundary, actuation :: AbstractIOSetup, input_signals :: AbstractArray{T}) where T <: Real  # in-place 
    α = λ/(c*ρ) # Diffusivity

    @inbounds for iz in 1:Nz
        @inbounds for iy in 1 : Ny
            @inbounds for ix in 2 : Nx-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                
                dθ[i] = α * (θ[i-1] + θ[i+1] - 2*θ[i])/Δx^2 
            end
            i1 = (iz-1)*Nx*Ny + (iy-1)*Nx + 1      # West
            i2 = (iz-1)*Nx*Ny + (iy-1)*Nx + Nx     # East
            
            emission_west = getEmission(boundary, i1, :west)
            emission_east = getEmission(boundary, i2, :east)
    
            char_west, in_idx_west = getActuation(actuation, i1, :west)
            char_east, in_idx_east = getActuation(actuation, i2, :east)

            ϕin_west = char_west * input_signals[in_idx_west]
            ϕin_east = char_east * input_signals[in_idx_east]

            dθ[i1] = α * (θ[i1+1] - θ[i1]) / Δx^2 + (emit(θ[i1], emission_west) + ϕin_west)/(c*ρ*Δx)
            dθ[i2] = α * (θ[i2-1] - θ[i2]) / Δx^2 + (emit(θ[i2], emission_east) + ϕin_east)/(c*ρ*Δx)
        end
    end
    nothing 
end



function diffusion_static_y!(dθ,θ,Nx, Ny, Nz, Δy, λ :: Real, c :: Real, ρ :: Real, boundary :: CubicBoundary, actuation :: AbstractIOSetup, input_signals :: AbstractArray{T}) where T <: Real # in-place
    α = λ/(c*ρ)
    
    @inbounds for iz in 1:Nz
        @inbounds for ix in 1 : Nx
            @inbounds for  iy in 2 : Ny-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                
                dθ[i] = dθ[i] + α * (θ[i-Nx] + θ[i+Nx] - 2*θ[i])/Δy^2 
    
            end
            i1 = (iz-1)*Nx*Ny + ix                 # South
            i2 = (iz-1)*Nx*Ny + Nx*(Ny-1) + ix     # North
            
            emission_south = getEmission(boundary, i1, :south)
            emission_north = getEmission(boundary, i2, :north)
    
            char_south, in_idx_south = getActuation(actuation, i1, :south)
            char_north, in_idx_north = getActuation(actuation, i2, :north)

            ϕin_south = char_south * input_signals[in_idx_south]
            ϕin_north = char_north * input_signals[in_idx_north]

            dθ[i1] = dθ[i1] + α * (θ[i1+Nx] - θ[i1]) / Δy^2 + (emit(θ[i1], emission_south) + ϕin_south)/(c*ρ*Δy)
            dθ[i2] = dθ[i2] + α * (θ[i2-Nx] - θ[i2]) / Δy^2 + (emit(θ[i2], emission_north) + ϕin_north)/(c*ρ*Δy)
        end
    end

    nothing 
end

function diffusion_static_z!(dθ,θ,Nx, Ny, Nz, Δz, λ :: Real, c :: Real, ρ :: Real, boundary :: CubicBoundary, actuation :: AbstractIOSetup, input_signals :: AbstractArray{T}) where T <: Real # in-place
    α = λ/(c*ρ)
    
    
    @inbounds for ix in 1 : Nx
        @inbounds for  iy in 1 : Ny
            @inbounds for iz in 2:Nz-1

                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                
                dθ[i] = dθ[i] + α * (θ[i-Nx*Ny] + θ[i+Nx*Ny] - 2*θ[i])/Δz^2 
    
            end
            
            i1 = (iy-1)*Nx + ix                 # Underside
            i2 = (Nz-1)*Nx*Ny + (iy-1)*Nx + ix  # Topside
            
            emission_underside = getEmission(boundary, i1, :underside)
            emission_topside   = getEmission(boundary, i2, :topside)
    
            char_underside, in_idx_underside = getActuation(actuation, i1, :underside)
            char_topside, in_idx_topside     = getActuation(actuation, i2, :topside)

            ϕin_underside = char_underside * input_signals[in_idx_underside]
            ϕin_topside   = char_topside   * input_signals[in_idx_topside]


            dθ[i1] = dθ[i1] + α * (θ[i1+Nx*Ny] - θ[i1]) / Δz^2 + (emit(θ[i1], emission_underside) + ϕin_underside)/(c*ρ*Δz)
            dθ[i2] = dθ[i2] + α * (θ[i2-Nx*Ny] - θ[i2]) / Δz^2 + (emit(θ[i2], emission_topside  ) + ϕin_topside  )/(c*ρ*Δz)
        end
    end

    nothing 
end





function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatrod :: HeatRod, 
                    property :: DynamicIsotropic, 
                    boundary :: CubicBoundary ,  
                    actuation :: AbstractIOSetup, 
                    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}
    
    λ(x) = specifyproperty(x, property.λ) # thermal conductivity
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    Nx = heatrod.heatcells
    Δx = heatrod.sampling
    
    diffusion_dynamic_x!(dθ, θ, Nx, 1, 1, Δx, λ, c, ρ, boundary, actuation, input_signals)
end


function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatplate :: HeatPlate, 
                    property :: DynamicIsotropic, 
                    boundary :: CubicBoundary,
                    actuation :: AbstractIOSetup, 
                    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}
   
    λ(x) = specifyproperty(x, property.λ) # thermal conductivity
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny) = heatplate.heatcells
    (Δx,Δy) = heatplate.sampling
        
    diffusion_dynamic_x!(dθ, θ, Nx, Ny, 1, Δx, λ, c, ρ, boundary, actuation, input_signals)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, 1, Δy, λ, c, ρ, boundary, actuation, input_signals)
end

function diffusion!(dθ :: AbstractArray{T1}, 
                    θ :: AbstractArray{T2},  
                    heatcuboid :: HeatCuboid, 
                    property :: DynamicIsotropic, 
                    boundary :: CubicBoundary,
                    actuation :: AbstractIOSetup, 
                    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}
   
    λ(x) = specifyproperty(x, property.λ) # thermal conductivity
    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling
        
    diffusion_dynamic_x!(dθ, θ, Nx, Ny, Nz, Δx, λ, c, ρ, boundary, actuation, input_signals)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, Nz, Δy, λ, c, ρ, boundary, actuation, input_signals)
    diffusion_dynamic_z!(dθ, θ, Nx, Ny, Nz, Δz, λ, c, ρ, boundary, actuation, input_signals)
end




function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatplate :: HeatPlate, 
    property :: DynamicAnisotropic, 
    boundary :: CubicBoundary,
    actuation :: AbstractIOSetup, 
    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}

    λx(x) = specifyproperty(x, property.λx) # thermal conductivity in x-direction
    λy(x) = specifyproperty(x, property.λy) # ... in y-direction

    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny) = heatplate.heatcells
    (Δx,Δy) = heatplate.sampling

    diffusion_dynamic_x!(dθ, θ, Nx, Ny, 1, Δx, λx, c, ρ, boundary, actuation, input_signals)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, 1, Δy, λy, c, ρ, boundary, actuation, input_signals)
end

function diffusion!(dθ :: AbstractArray{T1}, 
    θ :: AbstractArray{T2},  
    heatcuboid :: HeatCuboid, 
    property :: DynamicAnisotropic, 
    boundary :: CubicBoundary,
    actuation :: AbstractIOSetup, 
    input_signals :: AbstractArray{T3} ) where {T1 <: Real, T2 <: Real, T3 <: Real}
    
    λx(x) = specifyproperty(x, property.λx) # thermal conductivity in x-direction
    λy(x) = specifyproperty(x, property.λy) # ... in y-direction
    λz(x) = specifyproperty(x, property.λz) # ... in z-direction

    c(x) = specifyproperty(x, property.c) # specific heat capacity
    ρ(x) = specifyproperty(x, property.ρ) # mass density

    (Nx,Ny,Nz) = heatcuboid.heatcells 
    (Δx,Δy,Δz) = heatcuboid.sampling

    diffusion_dynamic_x!(dθ, θ, Nx, Ny, Nz, Δx, λx, c, ρ, boundary, actuation, input_signals)
    diffusion_dynamic_y!(dθ, θ, Nx, Ny, Nz, Δy, λy, c, ρ, boundary, actuation, input_signals)
    diffusion_dynamic_z!(dθ, θ, Nx, Ny, Nz, Δz, λz, c, ρ, boundary, actuation, input_signals)
end




function diffusion_dynamic_x!(dθ,θ,Nx, Ny, Nz, Δx, λ :: Function, c :: Function, ρ :: Function, boundary :: CubicBoundary, actuation :: AbstractIOSetup, input_signals :: AbstractArray{T}) where T <: Real 
    
    @inbounds for iz in 1:Nz
        @inbounds for iy in 1 : Ny
            @inbounds for ix in 2 : Nx-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                x̅1 =  (θ[i-1] + θ[i])/2
                x̅2 =  (θ[i+1] + θ[i])/2
                
                dθ[i] = (λ(x̅1) * θ[i-1] + λ(x̅2) * θ[i+1] - (λ(x̅1) + λ(x̅2))*θ[i])/(Δx^2 * ρ(θ[i]) * c(θ[i]))
            end
            i1 = (iz-1)*Nx*Ny + (iy-1)*Nx + 1      # West
            i2 = (iz-1)*Nx*Ny + (iy-1)*Nx + Nx     # East
         
            x̅1 =  (θ[i1] + θ[i1+1])/2
            x̅2 =  (θ[i2] + θ[i2-1])/2
            
            emission_west = getEmission(boundary, i1, :west)
            emission_east = getEmission(boundary, i2, :east)
    
            char_west, in_idx_west = getActuation(actuation, i1, :west)
            char_east, in_idx_east = getActuation(actuation, i2, :east)

            ϕin_west = char_west * input_signals[in_idx_west]
            ϕin_east = char_east * input_signals[in_idx_east]

            dθ[i1] = λ(x̅1) * (θ[i1+1] - θ[i1]) / (Δx^2 * ρ(θ[i1]) * c(θ[i1])) + (emit(θ[i1], emission_west) + ϕin_west)/(Δx * ρ(θ[i1]) * c(θ[i1]))
            dθ[i2] = λ(x̅2) * (θ[i2-1] - θ[i2]) / (Δx^2 * ρ(θ[i2]) * c(θ[i2])) + (emit(θ[i2], emission_east) + ϕin_east)/(Δx * ρ(θ[i2]) * c(θ[i2]))
        end
    end

    

    nothing 
end



function diffusion_dynamic_y!(dθ,θ,Nx, Ny, Nz, Δy, λ :: Function, c :: Function, ρ :: Function, boundary :: CubicBoundary, actuation :: AbstractIOSetup, input_signals :: AbstractArray{T}) where T <: Real 
    @inbounds for iz in 1:Nz
        @inbounds for ix in 1 : Nx
            @inbounds for  iy in 2 : Ny-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
    
    
                x̅1 =  (θ[i-Nx] + θ[i])/2
                x̅2 =  (θ[i+Nx] + θ[i])/2
                
                dθ[i] = dθ[i] +  (λ(x̅1) * θ[i-Nx] + λ(x̅2) * θ[i+Nx] - (λ(x̅1) + λ(x̅2))*θ[i])/(Δy^2 * ρ(θ[i]) * c(θ[i]))
    
            end
            i1 = (iz-1)*Nx*Ny + ix                 # South
            i2 = (iz-1)*Nx*Ny + Nx*(Ny-1) + ix     # North
    
            x̅1 =  (θ[i1] + θ[i1+Nx])/2
            x̅2 =  (θ[i2] + θ[i2-Nx])/2
    
            emission_south = getEmission(boundary, i1, :south)
            emission_north = getEmission(boundary, i2, :north)
    
            char_south, in_idx_south = getActuation(actuation, i1, :south)
            char_north, in_idx_north = getActuation(actuation, i2, :north)

            ϕin_south = char_south * input_signals[in_idx_south]
            ϕin_north = char_north * input_signals[in_idx_north]

            dθ[i1] = dθ[i1] + λ(x̅1) * (θ[i1+Nx] - θ[i1]) / (Δy^2 * ρ(θ[i1]) * c(θ[i1])) + (emit(θ[i1],emission_south) + ϕin_south)/ (Δy * ρ(θ[i1]) * c(θ[i1]))
            dθ[i2] = dθ[i2] + λ(x̅2) * (θ[i2-Nx] - θ[i2]) / (Δy^2 * ρ(θ[i2]) * c(θ[i2])) + (emit(θ[i2],emission_north) + ϕin_north)/ (Δy * ρ(θ[i2]) * c(θ[i2]))
        end
    end   

    nothing 
end




function diffusion_dynamic_z!(dθ,θ,Nx, Ny, Nz, Δz,  λ :: Function, c :: Function, ρ :: Function, boundary :: CubicBoundary, actuation :: AbstractIOSetup, input_signals :: AbstractArray{T}) where T <: Real 
    
    @inbounds for ix in 1 : Nx
        @inbounds for  iy in 1 : Ny
            @inbounds for iz in 2:Nz-1
                i = (iz-1)*Nx*Ny + (iy-1)*Nx + ix
                x̅1 =  (θ[i-Nx*Ny] + θ[i])/2
                x̅2 =  (θ[i+Nx*Ny] + θ[i])/2
                
                dθ[i] = dθ[i] +  (λ(x̅1) * θ[i-Nx*Ny] + λ(x̅2) * θ[i+Nx*Ny] - (λ(x̅1) + λ(x̅2))*θ[i])/(Δz^2 * ρ(θ[i]) * c(θ[i]))
            end
            
            i1 = (iy-1)*Nx + ix                 # Underside
            i2 = (Nz-1)*Nx*Ny + (iy-1)*Nx + ix  # Topside
            
            x̅1 =  (θ[i1] + θ[i1+Nx])/2
            x̅2 =  (θ[i2] + θ[i2-Nx])/2

            emission_underside = getEmission(boundary, i1, :underside)
            emission_topside   = getEmission(boundary, i2, :topside  )

            char_underside, in_idx_underside = getActuation(actuation, i1, :underside)
            char_topside, in_idx_topside     = getActuation(actuation, i2, :topside)

            ϕin_underside = char_underside * input_signals[in_idx_underside]
            ϕin_topside   = char_topside   * input_signals[in_idx_topside]

    
            dθ[i1] = dθ[i1] + λ(x̅1) * (θ[i1+Nx*Ny] - θ[i1]) / (Δz^2 * ρ(θ[i1]) * c(θ[i1])) + (emit(θ[i1],emission_underside) + ϕin_underside)/ (Δz * ρ(θ[i1]) * c(θ[i1]))
            dθ[i2] = dθ[i2] + λ(x̅2) * (θ[i2-Nx*Ny] - θ[i2]) / (Δz^2 * ρ(θ[i2]) * c(θ[i2])) + (emit(θ[i2],emission_topside)   + ϕin_topside  )/ (Δz * ρ(θ[i2]) * c(θ[i2]))
        end
    end

    nothing 
end
