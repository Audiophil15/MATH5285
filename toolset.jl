using Polynomials
using PNGFiles

function newton(p, z)
	return z-p(z)/(derivative(p)(z))
end

function convergence(polynomial, start, lim=limit, threshold=0.001)
	zn = [start]
	c = 0
	while c < lim && abs(zn[end]) < threshold
		append!(zn, newton(polynomial, zn[end]))
		c += 1
	end
	return c/lim
end

function convergencemap(polynomial, start, lim=1e1, threshold=0.001; verbose=false)
	zn = start
	zn2 = copy(zn)
	cr = zero(real.(zn))
	for c in 1:lim
		zn2 = copy(zn)
		if verbose
			println(c)
		end
		if verbose
			@time begin
				zn = newton.(polynomial, zn)
			end
		else
			zn = newton.(polynomial, zn)
		end
		cr += abs.(zn-zn2).>threshold
		if verbose
			println(zn[1,1])
			println()
		end
	end
	return (zn, cr)
end

function scalecr(cr, lim)
	return 1 .- sqrt.(cr/lim)
end

function juliasteps(z0, p, lim=30)
	t = [z0]
	for i in 1:lim
		push!(t, newton(p, t[end]))
	end
	return t
end

function zlisttocrcoord(t, ba, bb, steps)
	xi = real.(t)
	yi = imag.(t)
	xii = (xi .- ba) ./ (bb-ba) .* steps .+ 1
	yii = (yi .- ba) ./ (bb-ba) .* steps .+ 1
	return xii, yii
end

function bfsconnected(vmap, coord)
	if typeof(coord) == typeof(CartesianIndex(0,0))
		coord = collect(coord.I)
	end
	vv = []
	av = [coord]
	while length(av) != 0
		c = av[1]
		for i in [max(1, c[1]-1), min(vmap.size[1], c[1]+1)]
			for j in [max(1, c[2]-1), min(vmap.size[2], c[2]+1)]
				cc = [i,c[2]]
				if vmap[cc...] == 0. && !(cc in av) && !(cc in vv)
					push!(av, cc)
				end
				cc = [c[1],j]
				if vmap[cc...] == 0. && !(cc in av) && !(cc in vv)
					push!(av, cc)
				end
			end
		end
		push!(vv, popfirst!(av))
	end
	return vv
end

function savefigure(img, pathname="fig.png", savefig=true)
	if occursin("/", pathname)
		mkpath(join(split(pathname, "/")[begin:end-1], "/")*"/")
	end
	PNGFiles.save(pathname, img)
	return img
end

function createfigure(cr, pathname="fig.png", grayscale=false, savefig=true)
	img = 0
	if size(cr[1]) == ()
		global img = colorview(Gray, cr)
	else
		global img = colorview(RGB, cr)
	end

	if savefig
		savefigure(img, pathname, savefig)
	end
	return img
end


# function savecr(cr, pathname="cr.jld")
# end

function squaremap(xa, xb, ya, yb, stepsize; endpointx=true, endpointy=true)
	x = xa:stepsize:xb
	y = ya:stepsize:yb
	if !endpointx && xa != xb
		x = xa:stepsize:xb-stepsize
	end
	if !endpointy && ya != yb
		y = ya:stepsize:yb-stepsize
	end
	xy = [i+j*im for j in y, i in x]
	return xy
end
