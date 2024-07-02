using SymPy
using SpecialFunctions
using Latexify

@syms x, y, z
@syms a::(real, positive), b::(real, positive), c::(real, positive), n::(real, positive)
expr = x^n * exp(- a*x^2 - b*x - c)

# The actual Gaussian integral evaluation (From SymPy)

for i in 1:30
    currentExpr = expr.subs(n, i)
    currentExprInt = factor(integrate(currentExpr, (x, -oo, oo))).subs(erf(b/(2*sqrt(a))) + erfc(b/(2*sqrt(a))), 1)
    if i % 2 == 0
        display(currentExprInt.args[3])
    else
        display(currentExprInt.args[4])
    end
    display(currentExprInt)
end

# Testing part