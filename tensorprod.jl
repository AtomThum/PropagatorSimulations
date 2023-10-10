using Plots, Tensorial

simrange = 5 # Size of simulation (Because the way it's coded, ±simrange will have an infinite potential barrier)
dx = 0.02
simlength = 5
dt = 0.1
N = round(Int, simrange/(2 * dx)) # Number of x slots

# Position vector
x = [(n * dx) for n in -N:N] 

# range(start = -simrange, stop = simrange, length = N)
ψ = [exp(-(n * dx - 1)^2) for n = -N:N] # Gaussian wave packets no.1
φ = [exp(-(n * dx + 1)^2) for n = -N:N] # Gaussian wave packets no.2
ψφ = [exp(-(n * dx + 1)^2)*exp(-(n * dx - 1)^2) for n = -N:N]
# Probability density for plotting
probabilityDensity(Ψ) = real(Ψ)^2 + imag(Ψ)^2
plot(x, [ψ φ ψφ])