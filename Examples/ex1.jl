using Polynomials
using TypedPolynomials
using LinearAlgebra

p=Polynomials.Polynomial([1,2,3])
P=integrate(p)
P(1)-P(0)

@polyvar x y

monoms=TypedPolynomials.monomials((x, y), 0:2)
p=monoms[2]
p(x=>2,y=>3)
