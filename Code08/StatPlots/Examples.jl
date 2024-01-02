using StatsPlots, RDatasets  
mlmf = dataset("mlmRev","Gcsemv"); 

wF = collect(skipmissing(mlmf[mlmf.Gender .== "F", :Written]));  
wM = collect(skipmissing(mlmf[mlmf.Gender .== "M", :Written]));  
cF = collect(skipmissing(mlmf[mlmf.Gender .== "F", :Course]));  
cM = collect(skipmissing(mlmf[mlmf.Gender .== "M", :Course]));

#=
  By using the @df macro pass columns within an array and then call the density
  function to display the spectral density of the four reduced datasets.
=#

# Set up the legend
labs = ["Exam (Girls)" "Exam (Boys)" "Course (Girls)" "Course (Boys)"];  
 
@df mlmf density([wF, wM, cF, cM],labels=labs, legend = :topleft)

# Get data without any missing values in either of the two sets of marks

mlmfX = mlmf[completecases(mlmf),:] 

#=
  This reduces the original dataset from 1905 to 1523 values, still a reasonable
  number to work with. However, it could be argued that missing values were due
  to examinations missed and/or coursework not presented and so these values
  should be set to zero. 
=#
@df mlmfX plot(:School,[:Written :Course])
 
@df mlmfX corrplot([:Written :Course])
@df mlmfX dotplot([:Written :Course])

using Query 
mlmfX |> @filter(_.Written > 60) |> @df scatter(:Course)

