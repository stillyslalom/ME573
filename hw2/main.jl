using ProgressMeter # Load progress meter display library

## Note on performance:
# In Julia, if a non-constant variable is defined outside of a function and 
# then called by that function, performance decreases significantly, since the
# compiler must check the type of the variable for each iteration. Best 
# practice is setting all global vars as constants, then iterating within a 
# function so the compiler can optimize.


## Enter given dimensions =====================================================
const L_x 	= 8.		# Length of domain in x-direction
const L_y 	= 6.5		# Length of domain in y-direction
const x_c 	= 0.75*L_x;	# x-coord of 'pressure source'
const y_c 	= 0.35*L_y;	# y-coord of 'pressure source'
const A 	= 2.
const σ 	= 0.05*min(L_x,L_y) # Note: can use Greek characters (αβγδ, etc.)
const maxerr= 1e-8		# Maximum allowable solution error

## Generate grid ==============================================================
const nx 	= 80 			# Number of mesh points in the x-direction
const ny 	= 80 			# Number of mesh points in the y-direction
const sqnxny= sqrt(nx+ny)	# Precompute value (used in error measurement)
const dx 	= L_x/(nx-1); 	# Distance btwn adjacent mesh points in the x-dir
const dy 	= L_y/(ny-1);	# Distance btwn adjacent mesh points in the y-dir
# Precompute dx^2 and dy^2, since they are used frequently in calculations
const dx2 	= dx^2
const dy2 	= dy^2 	
# X[i], Y[j] is equivalent to x_i, y_j
X = linspace(0.0, L_x, nx) 	# Generate vector with X_i = dx*(i-1)
Y = linspace(0.0, L_y, ny) 	# Generate vector with Y_i = dy*(j-1)
XX = repmat(X , 1, ny) 	# Repeat the X-loc vector across all y-values
YY = repmat(Y', nx, 1)	# Repeat the Y-loc vector across all x-values

## Generate initial conditions ================================================
# Compute R as def'd by handout:
const R = [-A*exp(-((x-x_c)^2+(y-y_c)^2)/(2*σ)) for x = X, y = Y]
# Precomputes R-term for solution:
const Rmod = convert(Array{Float64,2},dx2*dy2*R/(2*(dx2+dy2)))		

const denom = 2*(dx2+dy2)	# Precomputes denominator of P-term for solution
P0 	 = zeros(nx,ny);		# Initialize pressure matrix as zero everywhere
for j = 2:ny-1				# Set P0 equal to R for all non-(P=0) boundaries
	for i = 2:nx
		P0[i,j] = Rmod[i,j]
	end
end
const P_init = copy(P0)		# Lock initial values of P into constant

## Method wrapper =============================================================
# Initializes variables and iterates solver until error is sufficiently small.
# Note that in Julia, multiple functions with different argument types can be 
# given the same name, and the compiler will choose the appropriate function
# at runtime

function methodwrap(solver, maxiter::Int64)
	P = copy(P_init)
	iter = 1 		# Set initial iteration count
	err  = 1.0 		# Initialize error value
	# Progress: displays progress bar at command prompt
	prog = Progress(maxiter,.05, "Solving using $solver: ", 10) 
	while (err > maxerr) && (iter < maxiter) # while(not converged)
		P, err = solver(P, maxiter)			 # iteratively solve for P
		next!(prog) 						 # increase progress bar counter
		iter += 1 							 # increase iteration count; end
	end
	println()
	return (P, iter, err)
end

function methodwrap(solver, maxiter::Int64, ω::Float64)
	P = copy(P_init)
	iter = 1 		# Set initial iteration count
	err  = 1.0 		# Initialize error value
	prog = Progress(maxiter,.05, "Solving using $solver with omega=$ω: ", 10)
	while (err > maxerr) && (iter < maxiter)
		P, err = solver(P, maxiter, ω)
		next!(prog) # Iterates progress bar counter
		iter += 1
	end
	println()
	return (P, iter, err)
end

## Function solvers ===========================================================
# The solvers themselves do not iterate--they are called by the iterating 
# function wrapper until convergence.
# Error is computed using the relative root mean square method: the Euclidean
# norm of the difference between successive solution matrices is itself 
# normalized to the square root of the sum of solution matrix dimensions.
# Descriptions of the solvers can be found in the report at:
# https://www.dropbox.com/s/x25q6eb9xu8otza/me573_hw2.pdf?dl=0

function jacobi(P::Array{Float64,2}, maxiter::Int64)
	P_old = copy(P)
	for j = 2:ny-1
		# Main body loop
		for i = 2:nx-1 
			 @inbounds P[i,j] = ((P_old[i+1,j] + P_old[i-1,j])*dx2 + (P_old[i,j+1] + P_old[i,j-1])*dy2)/denom-Rmod[i,j]
		end
		# RHS boundary condition loop (dP/dx=0)
		for i = nx
			P[i,j] = (2*P_old[i-1,j]*dx2 + (P_old[i,j+1] + P_old[i,j-1])*dy2)/denom-Rmod[i,j]
		end
	end
	err = vecnorm(P-P_old)/sqnxny 
	return (P, err)
end

function gauss_seidel(P::Array{Float64,2}, maxiter::Int64)
	P_old = copy(P)
	for j = 2:ny-1
		# Main body loop
		for i = 2:nx-1 
			 @inbounds P[i,j] = ((P[i+1,j] + P[i-1,j])*dx2 
			 				   + (P[i,j+1] + P[i,j-1])*dy2)/denom-Rmod[i,j]
		end
		# RHS boundary condition loop (dP/dx=0)
		for i = nx
			P[i,j] = (2*P[i-1,j]*dx2 + (P[i,j+1] + P[i,j-1])*dy2)/denom-Rmod[i,j]
		end
	end
	err = vecnorm(P-P_old)/sqnxny
	return (P, err)
end

function SOR(P::Array{Float64,2}, maxiter::Int64, ω::Float64)
	P_old = copy(P)
	for j = 2:ny-1
	# Main body loop
		for i = 2:nx-1 
			 @inbounds P[i,j] = ω*(((P[i+1,j] + P[i-1,j])*dx2 
			 				   + (P[i,j+1] + P[i,j-1])*dy2)/denom-Rmod[i,j]-P_old[i,j])+P_old[i,j]
		end
		# RHS boundary condition loop (dP/dx=0)
		for i = nx
			P[i,j] = ω*((2*P[i-1,j]*dx2 + (P[i,j+1] + P[i,j-1])*dy2)/denom-Rmod[i,j]-P_old[i,j])+P_old[i,j]
		end
	end
	err = vecnorm(P-P_old)/sqnxny
	return (P, err)
end

## Automatic method comparisons ===============================================
# These functions automatically run comparisons of solution methods and return 
# the results in a standard format (P_solution, iteration count, final error)

function method_compare(maxiter::Int64)
	methods = [jacobi, gauss_seidel]
	lm = length(methods)
	P = zeros(nx,ny,lm)
	it = zeros(Int64,lm)
	errval = zeros(lm)
	for i = 1:length(methods)
		@time P[:,:,i], it[i], errval[i] = methodwrap(methods[i], maxiter)
		println()
	end
	return P, it, errval
end

function method_compare(maxiter::Int64, omegas::Array{Float64,1})
	lo = length(omegas)
	P = zeros(nx,ny,lo)
	it = zeros(Int64,lo)
	errval = zeros(lo)
	for i in 1:lo
		@time P[:,:,i], it[i], errval[i] = methodwrap(SOR, maxiter, omegas[i])
		println()
	end
	return P, it, errval
end

function run_comparos()
	jacobi_gs_res = method_compare(30000);
	SOR_res		  = method_compare(30000,[1.3:.1:2.]);
	return (jacobi_gs_res, SOR_res)
end