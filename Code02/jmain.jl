# Program to create a Julia set on disk which can be run from command line.

# Include main function (jset.jl) and routine to create the disk file. 
# This does not require use of Images and ImgView packages and so is quicker.

include("jset.jl")
include("pgmfile.jl")

# Define image size and the 'c0' parameter

h = 400; w = 800; m = Array{Union{Nothing, Int}}(nothing, h, w);
c0 = -0.8+0.16im;

pgm_name = "julia.pgm";

# Time how long it will take to execute.

t0 = time();
for y=1:h, x=1:w
    c = complex((x-w/2)/(w/2), (y-h/2)/(w/2))
    m[y,x] = juliaset(c, c0, 256)
end
t1 = time();

# Save image on disk and output execution time.
   
create_pgmfile(m, pgm_name);
print("Written $pgm_name\nFinished in $(round((t1-t0),digits = 4)) seconds.\n");
