using Plots
using Polynomials
using Images, ImageView
using JLD

include("toolset.jl")

# p = Polynomial([-0.36,-0.64,0,1])
# p = Polynomial([-35028, 86868, -72419, -8509, 29746, -5842, -2363, 635])
# 635 x^7 - 2363 x^6 - 5842 x^5 + 29746 x^4 - 8509 x^3 - 72419 x^2 + 86868 x - 35028

limit = Int(1e2)
bay = -5
bby = 5
bax = -5
bbx = 5
steps = 1000
if steps%1000!=0
	println("There must be a multiple of 1000 for the steps.")
	exit(-1)
end

stepsize = 0.01 # (bby-bay)/steps

chunksx = ceil(Int, (bbx-bax)/stepsize/1000) + Int(bbx == bax)
chunksy = ceil(Int, (bby-bay)/stepsize/1000) + Int(bby == bay)

mua = 0
mub = 0
paramstepsize = 0.001

zn = 0
cr = 0
img = 0
xy = 0
npxy = 0
p = 0
cmap = 0
npmap = 0

for i in mua:paramstepsize:mub

	global p = Polynomial([-0.36,-0.64,0,1])
	poly = p+Polynomial([i])
	println(poly)

	# Small image version
	@time begin

		for sy in 0:chunksy-1
			for sx in 0:chunksx-1

				println("Doing chunk $(sy*chunksy+sx) of $(chunksx*chunksy-1)")

				ya = bay+((bby-bay)/chunksy)*sy
				yb = bay+((bby-bay)/chunksy)*(sy+1)
				xa = bax+((bbx-bax)/chunksx)*sx
				xb = bax+((bbx-bax)/chunksx)*(sx+1)

				# Grid to use
				global xy = squaremap(xa, xb, ya, yb, stepsize, endpointx=(sx==chunksx-1), endpointy=(sy==chunksy-1))
				global npxy = newton.(poly, xy)
				
				# Creating file name and path
				num = rpad("$i", length(split("$i", ".")[1])+1+length(split("$stepsize", ".")[1]), "0")
				name = "steps=$steps lim=$limit"
				if mua==mub
					path = "Figures/p=$(p.coeffs)/x=($bax,$bbx) y=($bay,$bby)/"
					name = "img "*name
				else
					path = "Anims/p=$(p.coeffs)/x=($bax,$bbx) y=($bay,$bby)/"
					name = "anim [mu=$num] "*name
				end
				if chunksx*chunksy > 1
					name *= " part[$sy, $sx]"
				end
				extpng = ".png"
				
				global cmap
				global npmap
				ma = map(x->(real(x)-xa)/(xb-xa+Int(xa==xb))+im*(imag(x)-ya)/(yb-ya+Int(ya==yb)), xy)
				ma = map(x->min(1, max(0, real(x)))+im*min(1, max(0, imag(x))), ma)
				mb = map(x->(real(x)-xa)/(xb-xa+Int(xa==xb))+im*(imag(x)-ya)/(yb-ya+Int(ya==yb)), npxy)
				mb = map(x->min(1, max(0, real(x)))+im*min(1, max(0, imag(x))), mb)
				cmap = map(x->RGB(real(x), 0, imag(x)), ma)
				savefigure(cmap, path*"cmap [full] "*name*extpng)
				cmap = map(x->RGB(real(x), 0, 0), ma)
				savefigure(cmap, path*"cmap [real] "*name*extpng)
				cmap = map(x->RGB(0, 0, imag(x)), ma)
				savefigure(cmap, path*"cmap [imag] "*name*extpng)
				npmap = map(x->RGB(real(x), 0, imag(x)), mb)
				savefigure(npmap, path*"npmap [full] "*name*extpng)
				npmap = map(x->RGB(real(x), 0, 0), mb)
				savefigure(npmap, path*"npmap [real] "*name*extpng)
				npmap = map(x->RGB(0, 0, imag(x)), mb)
				savefigure(npmap, path*"npmap [imag] "*name*extpng)

				# x0 = xy[juliabubbles[36][1]...]
				# t = juliasteps(x0, p)
				# xii, yii = zlisttocrcoord(t, ba, bb, steps)
				# c = 1.:length(t)
				# anim = @animate for i in 1:length(xii)
				# 	plot(x, y, img, size=(1500,1500))
				# 	plot!(xii[1:i], yii[1:i])
				# 	scatter!(xii[1:i], yii[1:i],markersize=10,markerstrokewidth=0.2,zcolor=c,color=cgrad(:bam, rev = true))
				# end
				# gif(anim, "gifs/julia steps $x0.gif", fps=3)

				# Creating convergence map
				@time begin
					global zn
					global cr
					zn, cr = convergencemap(poly, xy, limit)
				end
				global cr
				cr = scalecr(cr, limit) # Scaling

				# Saving picture
				global img
				img = createfigure(cr, path*name*extpng)
			end
		end

	end

end

div = map(x->collect(x.I), findall(x->x==0., cr))
juliabubbles = []
for d in div
	j1 = bfsconnected(cr, d)
	push!(juliabubbles, copy(j1))
	filter!(x->!(x in j1), div)
end
println(juliabubbles)

x0 = xy[juliabubbles[36][1]...]
t = juliasteps(x0, p)
xii, yii = zlisttocrcoord(t, ba, bb, steps)
c = 1.:length(t)
anim = @animate for i in 1:length(xii)
	plot(x, y, img, size=(1500,1500))
	plot!(xii[1:i], yii[1:i])
	scatter!(xii[1:i], yii[1:i],markersize=10,markerstrokewidth=0.2,zcolor=c,color=cgrad(:bam, rev = true))
end
gif(anim, "gifs/julia steps $x0.gif", fps=3)

for i in eachindex(div)
	println(i)
	ing = copy(img)
	j0 = bfsconnected(cr, div[i])
	for coord in j0
		ing[coord...] = RGB(i/length(div),1,0)
		imshow(ing)
		# readline()
	end
end

for i in eachindex(juliabubbles)
	println(i)
	j0 = bfsconnected(cr, juliabubbles[i][0])
	for coord in j0
		ing[coord...] = RGB(i/length(div),0,0.5)
	end
end

x = bax:stepsize:bbx
y = bay:stepsize:bby
for i in eachindex(juliabubbles)
	println(i)
	j1 = juliabubbles[i]
	out = plot(x, y, img, size=(750,750))
	for j in eachindex(j1)
		x0 = xy[j1[j]...]
		t = juliasteps(x0, p)
		c = 1.:length(t)
		# xii, yii = zlisttocrcoord(t, bax, bbx, steps)
		# out = plot!(xii, yii, primary=false, linecolor=cgrad(:matter, length(j1), categorical = true)[j])
		# out = scatter!(xii, yii, markersize=10,markerstrokewidth=0.2,zcolor=c,color=cgrad(:bam, rev = true), primary=false)
		out = plot!(real.(t), imag.(t), primary=false, linecolor=cgrad(:matter, length(j1), categorical = true)[j])
		out = scatter!(real.(t), imag.(t), markersize=5,markerstrokewidth=0.2,zcolor=c,color=cgrad(:bam, rev = true), primary=false)
		# c = 1.:length(t)
		# anim = @animate for i in 1:length(xii)
		# 	plot(x, y, img, size=(1500,1500))
		# 	plot!(xii[1:i], yii[1:i])
		# 	scatter!(xii[1:i], yii[1:i],markersize=10,markerstrokewidth=0.2,zcolor=c,color=cgrad(:bam, rev = true))
		# end
		# gif(anim, "gifs/julia steps $x0.gif", fps=3)
	end
	savefig(out, "plots/path groups - group $i.png")
end

# 7, 26, 156, 13, 181, 210


