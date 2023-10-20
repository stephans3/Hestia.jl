
function diffusion_fast!(
    dθ::AbstractArray{T1},
    θ::AbstractArray{T2},
    heatplate::HeatPlate,
    property::StaticIsoProperty,
    boundary::CubicBoundary,
    actuation::AbstractIOSetup,
    input_signals::AbstractArray{T3},
) where {T1<:Real,T2<:Real,T3<:Real}

    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx, Ny) = heatplate.heatcells
    (Δx, Δy) = heatplate.sampling

    diffusion_static_x_fast!(
        dθ,
        θ,
        Nx,
        Ny,
        1,
        Δx,
        λ,
        c,
        ρ,
        boundary,
        actuation,
        input_signals,
    )
    diffusion_static_y_fast!(
        dθ,
        θ,
        Nx,
        Ny,
        1,
        Δy,
        λ,
        c,
        ρ,
        boundary,
        actuation,
        input_signals,
    )
end


function diffusion_static_x_fast!(
    dx,
    x,
    Nx,
    Ny,
    Nz,
    Δx,
    λ::Real,
    c::Real,
    ρ::Real,
    boundary::CubicBoundary,
    actuation::AbstractIOSetup,
    input_signals::AbstractArray{T},
) where {T<:Real}  # in-place 
    α = λ / (c * ρ) # Diffusivity

    @inbounds for iz = 1:Nz
        @inbounds for iy = 1:Ny
            @inbounds for ix = 2:Nx-1
                i = (iz - 1) * Nx * Ny + (iy - 1) * Nx + ix

                dx[i] = α * (x[i-1] + x[i+1] - 2 * x[i]) / Δx^2
            end
            i1 = (iz - 1) * Nx * Ny + (iy - 1) * Nx + 1      # West
            i2 = (iz - 1) * Nx * Ny + (iy - 1) * Nx + Nx     # East

            emission_west = getEmission(boundary, i1, :west)
            emission_east = getEmission(boundary, i2, :east)

            char_west, in_idx_west = getActuation(actuation, i1, :west)
            char_east, in_idx_east = getActuation(actuation, i2, :east)

            ϕin_west = char_west * input_signals[in_idx_west]
            ϕin_east = char_east * input_signals[in_idx_east]

            dx[i1] =
                α * (x[i1+1] - x[i1]) / Δx^2 +
                (emit(x[i1], emission_west) + ϕin_west) / (c * ρ * Δx)
            dx[i2] =
                α * (x[i2-1] - x[i2]) / Δx^2 +
                (emit(x[i2], emission_east) + ϕin_east) / (c * ρ * Δx)
        end
    end
    nothing
end



function diffusion_static_y_fast!(
    dx,
    x,
    Nx,
    Ny,
    Nz,
    Δy,
    λ::Real,
    c::Real,
    ρ::Real,
    boundary::CubicBoundary,
    actuation::AbstractIOSetup,
    input_signals::AbstractArray{T},
) where {T<:Real} # in-place
    α = λ / (c * ρ)

    @inbounds for iz = 1:Nz
        @inbounds for ix = 1:Nx
            @inbounds for iy = 2:Ny-1
                i = (iz - 1) * Nx * Ny + (iy - 1) * Nx + ix

                dx[i] = dx[i] + α * (x[i-Nx] + x[i+Nx] - 2 * x[i]) / Δy^2

            end
            i1 = (iz - 1) * Nx * Ny + ix                 # South
            i2 = (iz - 1) * Nx * Ny + Nx * (Ny - 1) + ix     # North

            emission_south = getEmission(boundary, i1, :south)
            emission_north = getEmission(boundary, i2, :north)

            char_south, in_idx_south = getActuation(actuation, i1, :south)
            char_north, in_idx_north = getActuation(actuation, i2, :north)

            ϕin_south = char_south * input_signals[in_idx_south]
            ϕin_north = char_north * input_signals[in_idx_north]

            dx[i1] =
                dx[i1] +
                α * (x[i1+Nx] - x[i1]) / Δy^2 +
                (emit(x[i1], emission_south) + ϕin_south) / (c * ρ * Δy)
            dx[i2] =
                dx[i2] +
                α * (x[i2-Nx] - x[i2]) / Δy^2 +
                (emit(x[i2], emission_north) + ϕin_north) / (c * ρ * Δy)
        end
    end

    nothing
end

function diffusion_fast_xy!(
    dθ::AbstractArray{T1},
    θ::AbstractArray{T2},
    heatplate::HeatPlate,
    property::StaticIsoProperty,
    boundary::CubicBoundary,
    actuation::AbstractIOSetup,
    input_signals::AbstractArray{T3},
) where {T1<:Real,T2<:Real,T3<:Real}

    λ = property.λ # Thermal conductivity
    c = property.c # capacity
    ρ = property.ρ # density    

    (Nx, Ny) = heatplate.heatcells
    (Δx, Δy) = heatplate.sampling

    diffusion_static_xy_fast!(
        dθ,
        θ,
        Nx,
        Ny,
        Δx,
        Δy,
        λ,
        c,
        ρ,
        boundary,
        actuation,
        input_signals,
    )
end


function diffusion_static_xy_fast!(
    dx,
    x,
    Nx,
    Ny,
    Δx,
    Δy,
    λ::Real,
    c::Real,
    ρ::Real,
    boundary::CubicBoundary,
    actuation::AbstractIOSetup,
    input_signals::AbstractArray{T},
) where {T<:Real}  # in-place 
    α = λ / (c * ρ)
    @inbounds for ix = 2:Nx-1, iy = 2:Ny-1
        i = (iy - 1) * Nx + ix
        dx[i] =
            α *
            ((x[i-1] + x[i+1] - 2 * x[i]) / Δx^2 + (x[i-Nx] + x[i+Nx] - 2 * x[i]) / Δy^2)
    end


    # West ix = 1
    @inbounds for iy = 2:Ny-1
        i = (iy - 1) * Nx + 1
        em = getEmission(boundary, i, :west)
        ind_ch, ind_idx = getActuation(actuation, i, :west)
        ϕin = ind_ch * input_signals[ind_idx]
        dx[i] =
            α * (x[i+1] - x[i]) / Δx^2 +
            (emit(x[i], em) + ϕin) / (c * ρ * Δx) +
            α * (x[i-Nx] + x[i+Nx] - 2 * x[i]) / Δy^2
    end

    # East ix = Nx
    @inbounds for iy = 2:Ny-1
        i = (iy - 1) * Nx + Nx
        em = getEmission(boundary, i, :east)
        ind_ch, ind_idx = getActuation(actuation, i, :east)
        ϕin = ind_ch * input_signals[ind_idx]
        dx[i] =
            α * (x[i-1] - x[i]) / Δx^2 +
            (emit(x[i], em) + ϕin) / (c * ρ * Δx) +
            α * (x[i-Nx] + x[i+Nx] - 2 * x[i]) / Δy^2
    end

    # South: iy = 1
    @inbounds for ix = 2:Nx-1
        i = ix
        em = getEmission(boundary, i, :south)
        ind_ch, ind_idx = getActuation(actuation, i, :south)
        ϕin = ind_ch * input_signals[ind_idx]
        dx[i] =
            α * (x[i+Nx] - x[i]) / Δy^2 +
            (emit(x[i], em) + ϕin) / (c * ρ * Δy) +
            α * (x[i-1] + x[i+1] - 2 * x[i]) / Δx^2
    end

    # North: iy = Ny
    @inbounds for ix = 2:Nx-1
        i = (Ny - 1) * Nx + ix
        em = getEmission(boundary, i, :north)
        ind_ch, ind_idx = getActuation(actuation, i, :north)
        ϕin = ind_ch * input_signals[ind_idx]
        dx[i] =
            α * (x[i-Nx] - x[i]) / Δy^2 +
            (emit(x[i], em) + ϕin) / (c * ρ * Δy) +
            α * (x[i-1] + x[i+1] - 2 * x[i]) / Δx^2
    end

    @inbounds begin
        # South West
        i = 1
        em_1 = getEmission(boundary, i, :west)
        ind_ch_1, ind_idx_1 = getActuation(actuation, i, :west)
        ϕin_1 = ind_ch_1 * input_signals[ind_idx_1]

        em_2 = getEmission(boundary, i, :south)
        ind_ch_2, ind_idx_2 = getActuation(actuation, i, :south)
        ϕin_2 = ind_ch_2 * input_signals[ind_idx_2]

        dx[i] =
            α * (x[i+1] - x[i]) / Δx^2 +
            (emit(x[i], em_1) + ϕin_1) / (c * ρ * Δx) +
            α * (x[i+Nx] - x[i]) / Δy^2 +
            (emit(x[i], em_2) + ϕin_2) / (c * ρ * Δy)

        # South East
        i = Nx
        em_1 = getEmission(boundary, i, :east)
        ind_ch_1, ind_idx_1 = getActuation(actuation, i, :east)
        ϕin_1 = ind_ch_1 * input_signals[ind_idx_1]

        em_2 = getEmission(boundary, i, :south)
        ind_ch_2, ind_idx_2 = getActuation(actuation, i, :south)
        ϕin_2 = ind_ch_2 * input_signals[ind_idx_2]

        dx[i] =
            α * (x[i-1] - x[i]) / Δx^2 +
            (emit(x[i], em_1) + ϕin_1) / (c * ρ * Δx) +
            α * (x[i+Nx] - x[i]) / Δy^2 +
            (emit(x[i], em_2) + ϕin_2) / (c * ρ * Δy)


        # North West
        i = (Ny - 1) * Nx + 1
        em_1 = getEmission(boundary, i, :west)
        ind_ch_1, ind_idx_1 = getActuation(actuation, i, :west)
        ϕin_1 = ind_ch_1 * input_signals[ind_idx_1]

        em_2 = getEmission(boundary, i, :north)
        ind_ch_2, ind_idx_2 = getActuation(actuation, i, :north)
        ϕin_2 = ind_ch_2 * input_signals[ind_idx_2]

        dx[i] =
            α * (x[i+1] - x[i]) / Δx^2 +
            (emit(x[i], em_1) + ϕin_1) / (c * ρ * Δx) +
            α * (x[i-Nx] - x[i]) / Δy^2 +
            (emit(x[i], em_2) + ϕin_2) / (c * ρ * Δy)

        # North East
        i = Ny * Nx
        em_1 = getEmission(boundary, i, :east)
        ind_ch_1, ind_idx_1 = getActuation(actuation, i, :east)
        ϕin_1 = ind_ch_1 * input_signals[ind_idx_1]

        em_2 = getEmission(boundary, i, :north)
        ind_ch_2, ind_idx_2 = getActuation(actuation, i, :north)
        ϕin_2 = ind_ch_2 * input_signals[ind_idx_2]

        dx[i] =
            α * (x[i-1] - x[i]) / Δx^2 +
            (emit(x[i], em_1) + ϕin_1) / (c * ρ * Δx) +
            α * (x[i-Nx] - x[i]) / Δy^2 +
            (emit(x[i], em_2) + ϕin_2) / (c * ρ * Δy)
    end

end
