using Plots
using Polynomials

const limit = 10e2

function newton(p, z)
	return z-p(z)/(derivative(p)(z))
end

p = Polynomial([0,2,3])

z = [4.0+2im]
c = 0
while c < limit && abs(z[end]) > 0.001
	append!(z, newton(p, z[end]))
end

println(z)

plot()

if isa(z[1],Complex)
	for v in z
		scatter!([v.re], [v.im])
	end
else
	x = -5:0.01:5
	plot!(x, p.(x))
	scatter!(z, p.(z))
end

gui()
readline()
