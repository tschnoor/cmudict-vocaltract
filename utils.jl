"""
Computes the frequency response, imput impedance,
and radiation impedance using a transmission line
model of the vocal tract. The model is a recreation
of Dr. Brad Story's MATLAB code which performs the
same function. His model is based on that described in 
[Sondhi and Schroeter (1987)](https://doi.org/10.1109/TASSP.1987.1165240).

The arguments to the function are:
- `a`: the area function which approximates the vocal tract (area as a function of distance from glottis)
- `l`: the length of each tubelet in the area function
- `cutoff`: the hard walled tube to extend the vocal tract (only used for impedance studies and normally set to number of tubelets+1)
- `max_f`: the ceiling of the frequency response
- `r`: the ratio of yielding wall resistance to mass
- `fw`: mechanical resonance frequency of the wall
- `ft`: lowest resonant frequency of the tract when closed on both ends 
- `q`: correction for thermal conductivity and viscosity normally set to 4 rad/s
- `c`: speed of sound
- `ρ`: density of air
"""
function transmission_line_model(
    a::Vector{Float64}=fill(1.0,44),
    l::Vector{Float64}=fill(.4,44),
    cutoff::Int64=length(a)+1,
    max_f::Int64=5000,
    r::Int64=408,
    fw::Int64=15,
    ft::Int64=200,
    q::Float64=4.0,
    c::Float64=35000.0,
    ρ::Float64=0.00114,
)

    # Preliminary calculations
    n_sections = length(a) # number of sections in area function
    n_points = Int64(round(max_f/5)) # maximum points in frequency spectrum
    delta_f = max_f/n_points
    frequencies = [delta_f:delta_f:delta_f*n_points;] # the frequencies at which response and impedance are calculated
    n_frequencies = length(frequencies) # number of frequencies at which ...
    ω = (2*π) .* frequencies # TODO define the rest below
    α = zeros(n_frequencies, 1) + im*(q*ω) 
    α = sqrt.(α)
    temp1 = r*ones(n_frequencies, 1) + im*(ω)
    temp2 = zeros(n_frequencies, 1) + im*(ω)
    den = temp1 .* temp2
    temp1 = ((2*π*fw)^2)*ones(n_frequencies, 1) + im*zeros(n_frequencies, 1)
    den = den + temp1
    num = zeros(n_frequencies, 1) + im*ω*(2*π*ft)^2
    β = (num ./ den) + α

    # print("Preliminary calculations complete")

    # Main loop
    # These are the ABCD matrices that make up the chain matrix for the vocal tract.
    # They are based on electrical transmission line theory.

    A = ones(n_points, 1);
    C = zeros(n_points, 1);
    B = zeros(n_points, 1);
    D = ones(n_points, 1);

    for k in 1:n_sections
    # for k in 1:3
        # println("iteration ", k)

        temp1 = r*ones(n_frequencies, 1) + im*(ω)
        temp2 = β + (zeros(n_frequencies, 1) + im*ω)
        γ = sqrt.(temp1 ./ temp2)
        σ = γ .* temp2

        nA = cosh.((l[k]/c) * σ)
        temp1 = -ρ*c/a[k] * γ
        nB = temp1 .* sinh.( (l[k]/c) * σ)
        temp1 = ones(length(γ), 1) ./ γ
        temp2 = -a[k]/(ρ*c) * temp1
        nC = temp2 .* sinh.((l[k]/c) * σ)
        nD = nA

        M1 = [A B; C D]
        M2 = [nA nB; nC nD]

        # println("M1: ", M1)
        # println(size(M1))
        # println("M2: ", M2)
        # println(size(M2))

        # M1 = vcat([A B], [C D])
        # M2 = vcat([nA nB], [nC nD])

        A = M2[1:n_points, 1].*M1[1:n_points, 1] + M2[1:n_points, 2].*M1[n_points+1:end, 1]
        B = M2[1:n_points, 1].*M1[1:n_points, 2] + M2[1:n_points, 2].*M1[n_points+1:end, 2]
        C = M2[n_points+1:end, 1].*M1[1:n_points, 1] + M2[n_points+1:end, 2].*M1[n_points+1:end, 1]
        D = M2[n_points+1:end, 1].*M1[1:n_points, 2] + M2[n_points+1:end, 2].*M1[n_points+1:end, 2]

        # println("A: ",A)
        # println("B: ",B)
        
    end
    # println("Finished loop")

    # Final calculations
    R = 128 .* ρ*c/(9*π^2 * a[n_sections])    
    L = 8 .* ρ*c/((3*π*c)*sqrt.(a[n_sections]*π)) 	
    
    temp1 = zeros(n_points, 1) + im*R*L*ω
    temp2 = R*ones(n_points, 1) + im*L*ω
    Zrad = temp1 ./ temp2

    temp1 = (Zrad .* D) - B
    temp2 = A - (C .* Zrad)
    z = temp1 ./ temp2
    h = Zrad ./ temp2
    f = frequencies

    S = temp2 ./ temp1

    h_alt = (-C .* B) ./ A+D

    # println("Calculations complete")

    return f, h, z, Zrad

end

struct Sequence
    text::String
    phones::Vector{Any}
    areas::Vector{Any}
    frequencies::Vector{Float64}
    responses::Vector{Any}
end

function create_structure(
    text::String,
    dict_file::String="cmudict.dict",
    area_file::String="areas.dict",
    tubelet_length::Float64=0.396825,
)

    doc = readlines(dict_file)
    entries = split.(doc, " ")
    phone_dict = Dict()
    for i in 1:length(entries)
        phone_dict[entries[i][1]] = entries[i][2:end]
    end

    doc = readlines(area_file)
    entries = split.(doc, " ")
    area_dict = Dict()
    for i in 1:length(entries)
        if entries[i][2:end] != ["NA"]
            af = parse.(Float64, entries[i][2:end])
            area_dict[entries[i][1]] = af
        end
    end

    words = split(text, " ")
    words = lowercase.(words)
    words = replace.(words, r"[.,;:!?]" => s"")
    phones = []
    for i in 1:length(words)
        push!(phones, phone_dict[words[i]])
    end

    areas = []
    for i in 1:length(phones)
        temp = []
        for j in 1:length(phones[i])
            key = replace(phones[i][j], r"[0-9]" => s"")
            if haskey(area_dict, key)
                push!(temp, area_dict[key])
            else
                push!(temp, nothing)
            end
        end
        push!(areas, temp)
    end

    frequencies = []
    responses = []
    for i in 1:length(areas)
        temp = []
        for j in 1:length(areas[i])
            if areas[i][j] != nothing
                af = areas[i][j]
                frequencies, h, _, _ = transmission_line_model(af, fill(tubelet_length, length(af)))
                push!(temp, 20 .* log10.(abs.(h)))
            else
                push!(temp, nothing)
            end
        end
        push!(responses, temp)
    end

    return Sequence(text, phones, areas, frequencies, responses)
end

function plot_vt_radii(
    areas::Vector{Float64},
    lengths::Vector{Float64}=fill(0.396825, length(areas)),
)

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
    a = @. sqrt(areas / π)
    l = lengths
    p = plot(xlabel="Length (cm)", ylabel="Tubelet Radius (cm)")
    x_offset = 0.0
    half_max_a = maximum(a) / 2
    for i in 1:length(a)
        y_offset = half_max_a - (a[i] / 2)
        plot!(rectangle(l[i], a[i], x_offset, y_offset), label="", color=:gray)
        x_offset += l[i]
    end
    return p
end