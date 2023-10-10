using LinearAlgebra
using MultivariatePolynomials, TypedPolynomials
using Plots

@polyvar x

# Propagator for the free particle (Let ℏ = 1)
K(x₀, x, t, m) = sqrt(m/(2*π*im*t)) * exp(-m*(x₀ - x)^2/t) # Scalar number
P(Ψ) = real.(Ψ)^2 + imag.(Ψ)^2 # Convert to probability distribution

simrange = 5 # Size of simulation (Because the way it's coded, ±simrange will have an infinite potential barrier)
dx = 0.01
simlength = 5
dt = 0.1
N = round(Int, simrange/(2 * dx)) # Number of x slots

posVect = [(n * dx) for n in -N:N] # Position basis
initVect = [exp(-(n * dx)^2 + 1) for n = -N:N] # Initial ket
plot(posVect, initVect) # Plotting the initial condition

# Propagator sampling
function propagate(initVect, targetTime, posVect)
    vectSize = size(posVect, 1)
    propVect = zeros(vectSize)im
    for i in range(1, vectSize)
        pdf = 0.0 + 0.0*im
        for j in range(1, vectSize)
            pdf = pdf + initVect[j]*K(posVect[j], posVect[i], targetTime, 1)
        end
        propVect[i] = pdf
    end
    return(propVect)
end

tValues = range(0, simlength, length = round(Int, simlength/dt))
animation = @animate for t in tValues
    global posVect
    global initVect
    newVect = normalize(P.(propagate(initVect, t, posVect)))
    plot(posVect, [newVect], title = "t = $t", yrange = (0, 0.25))
end
gif(animation, "freeparticle.gif", fps = 30)