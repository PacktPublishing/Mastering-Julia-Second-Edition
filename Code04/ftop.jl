#! /Users/malcolm/bin/julia -q
#
nargs = size(ARGS)[1]
if nargs > 2
  println("usage: ftop [dir [nfiles]]")
  exit(-1)
else
  dir = (nargs > 0) ? ARGS[1] : pwd()
  nf  = (nargs > 1) ? ARGS[2] : 20
end
efind =  `find $dir -type f -iname "*" -exec du -sh "{}" + `
println("Directory: $dir")
try
  run(pipeline(efind,`sort -rh`,`head -$nf`))
  println("Done.")
catch
  println("No files found.")
end
