using LinearAlgebra
using MultivariatePolynomials, TypedPolynomials
using Plots

@polyvar x

# Propagator for the free particle (Let ℏ = 1)
K(x₀, x, t, m) = sqrt(m/(2*π*im*t)) * exp(-m*(x₀ - x)^2/t) # Scalar number

simrange = 5 # Size of simulation (Because the way it's coded, ±simrange will have an infinite potential barrier)
dx = 0.02
simlength = 5
dt = 0.1
N = round(Int, simrange/(2 * dx)) # Number of x slots

# Initial condition
posVect = [-(n * dx) for n = -N:N]
# range(start = -simrange, stop = simrange, length = N)
initVect = [exp(-(n * dx)^2) for n = -N:N] # Gaussian wave packets

# Probability density for plotting
probabilityDensity(Ψ) = real(Ψ)^2 + imag(Ψ)^2

# Single point propagator
function propagatePoint(initVect, posVect, targetPos, targetTime, mass)
    pdf = 0
    vectSize = size(posVect, 1)
    for i in range(1, vectSize)
        initProb = initVect[i]
        initPos = posVect[i]
        pdf =+ initProb * K(initPos, targetPos, targetTime, mass)
    end
    return probabilityDensity(pdf)
end

# Propagating the whole wave function through time
function propagateVect(initVect, posVect, targetTime, mass)
    vectSize = size(posVect, 1)
    propagatedVect = zeros(vectSize)
    for j in range(1, vectSize)
        propagatedVect[j] = propagatePoint(
            initVect, posVect, posVect[j], targetTime, mass
        )
    end
    return propagatedVect
end

tValues = range(0, simlength, length = round(Int, simlength/dt))
animation = @animate for t in tValues
    global posVect
    global initVect
    newVect = normalize(propagateVect(initVect, posVect, t, 1))
    plot(posVect, [newVect], title = "t = $t")
end

gif(animation, "freeparticle.gif", fps = 15)