using Plots

simrange = 5 # Size of simulation (No potential barrier at edge), it just completely doesn't render at all.
dx = 0.02
simlength = 5
dt = 0.1
N = round(Int, simrange/(2 * dx)) # Number of x slots

# Propagator for single particle
K(x₀, x, t, m) = sqrt(m/(2*π*im*t)) * exp(-m*(x₀ - x)^2/t)
P(Ψ) = real.(Ψ)^2 + imag.(Ψ)^2

# Position vector
posVect = [(n * dx) for n in -N:N]

# range(start = -simrange, stop = simrange, length = N)
initVect1 = [exp(-(n * dx - 1)^2) for n = -N:N] # Gaussian wave packets no.1
initVect2 = [exp(-(n * dx + 1)^2) for n = -N:N] # Gaussian wave packets no.2

# Probability density for plotting
probabilityDensity(Ψ) = real(Ψ)^2 + imag(Ψ)^2
tensorProdPosition(Ψ₁, Ψ₂) = Ψ₁ .* Ψ₂
initVectCombined = tensorProdPosition(initVect1, initVect2)

Kn(targetX, x₁, x₂, t, m) = K(x₁, targetX, t, m)*K(x₂, targetX, t, m)

function doublePropagate(initVect, partOne, partTwo, targetTime, posVect)
    vectSize = size(posVect, 1)
    propVect = zeros(vectSize)im
    for i in range(1, vectSize)
        pdf = 0
        for x₁ in range(1, vectSize)
            for x₂ in range(1, vectSize)
                pdf += partOne[x₁]*partTwo[x₂]*Kn(posVect[i], posVect[x₁], posVect[x₂], targetTime, 1)
            end
        end
        propVect[i] = pdf
    end
    return propVect
end

tValues = range(0, simlength, length = round(Int, simlength/dt))
animation = @animate for t in tValues
    global posVect
    global initVect1
    global initVect2

    propVect = normalize(P.(doublePropagate(initVectCombined, initVect1, initVect2, t, posVect)))
    plot(posVect, [propVect], title = "t = $t", yrange = (0, 0.25))
end
gif(animation, "twofreeparticles.gif", fps = 30)