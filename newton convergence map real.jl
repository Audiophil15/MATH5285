using Plots
using Polynomials
using Images, ImageView
using PNGFiles

include("toolset.jl")

include("constants.jl")

p = Polynomial([-0.36,-0.64,0,1])
println("Working on $p")

limit = 1e7
ba = -4
bb = 0
steps = 1000
step = (bb-ba)/steps
chunks = Int(steps/100)

# img = 0
savefigs = true

# if steps < 10000
	# Small image version
	@time begin
			x = ba:step:bb
			
			@time begin
				xn, cr = convergencemap(p, x, limit)
			end

			cr = scalecr(cr, limit)

			img = plot(x, cr)
			if savefigs
				savefig(img, "plots/plot ($ba, $bb) $steps.svg")
			end
	end
gui()

# else
# 	# Big image version
# 	@time begin
# 		for s in 0:chunks-1
# 			baa = ba+((bb-ba)/chunks)*s
# 			bbb = ba+((bb-ba)/chunks)*(s+1)
# 			ba = baa
# 			bb = bbb
# 			x = (ba:step:bb)[begin:end-1]
# 			y = ba:step:bb
# 			xy = [i+j*im for j in y, i in x]
# 			println("Doing chunk $s")

# 			@time begin
# 				xn, cr = convergencemap(p, xy, limit)
# 			end

# 			cr = scalecr(cr, limit)

# 			img = colorview(Gray, cr)
# 			PNGFiles.save("img - $ba:$bb - $steps.png", img)
# 		end
# 	end

# end

# div = findall(x->x==0., cr)
# x0 = xy[div[4208]]
# t = juliasteps(x0)
# xii, yii = zlisttocrcoord(t, ba, bb, steps)
# c = 1.:length(t)
# anim = @animate for i in 1:length(xii)
# 	plot(img, size=(1500,1500))
# 	plot!(xii[1:i], yii[1:i])
# 	scatter!(xii[1:i], yii[1:i],markersize=10,markerstrokewidth=0.2,zcolor=c,color=cgrad(:bam, rev = true))
# end
# gif(anim, "gifs/julia steps $x0.gif", fps=3)

# imshow(img)

