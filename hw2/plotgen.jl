##
using PyPlot
using LaTeXStrings
include("main.jl")
pressure = methodwrap(SOR, 10000, 1.7)[1] # Get pressure solution array
#display_P(run_comparos());
function solnp()
	solnplot = contourf(XX, YY, pressure, 10, cmap = ColorMap("gray"))
	cbar = colorbar(solnplot)
	xlabel(L"x")
	ylabel(L"y")
	savefig("plots/solution_plot.svg")
end

function rp()
	rplot = contourf(XX, YY, R, 10, cmap = ColorMap("gray"))
	cbar2 = colorbar(rplot)
	xlabel(L"x")
	ylabel(L"y")
	savefig("plots/R_plot.svg")
end

function display_P(sol_tuples) # Fancy output of solution parameters
	jacobi_gs 	= sol_tuples[1]
	SOR_res 	= sol_tuples[2]
	is 		= [60 20 70]
	js 		= [16 40 60]
	for i = 1:2
		P 		= jacobi_gs[1][:,:,i]
		iter 	= jacobi_gs[2][i]
		errval 	= jacobi_gs[3][i]
		solver 	= jacobi_gs[4][i]
		println("Using the $solver solver, an error value of $errval")
		println("was reached after $iter iterations.")
		for j = 1:3
			println(" At i = $(is[j]), j = $(js[j]):  P = $(P[is[j],js[j]])")
		end
		println()
	end

	for i = 1:length(SOR_res[2])
		P 		= SOR_res[1][:,:,i]
		iter 	= SOR_res[2][i]
		errval 	= SOR_res[3][i]
		omega 	= SOR_res[4][i]
		println("Using the SOR solver with omega = $omega, an error value of")
		println("$errval was reached after $iter iterations.")
		for j = 1:3
			println(" At i = $(is[j]), j = $(js[j]):  P = $(P[is[j],js[j]])")
		end
		println()
	end
	return sol_tuples
end