function ftop(dir=".";nf=10)
  @assert nf > 0
  efind =  `find $dir -type f -iname "*" -exec du -sh "{}" + `
  try
    run(pipeline(efind,`sort -rh`,`head -$nf`))
    println("Done.")
  catch
    println("No files found.")
  end
end
