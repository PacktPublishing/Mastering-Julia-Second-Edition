## Mastering Julia - Second Edition

This is the code repository for Mastering Julia - Second Edition, published by Packt.


### What is this book about?

The first edition of this book was written when Julia was then relatively new (v0.4.0) and many changes and additions have been made in the subsequent years and this edition brings it up to date (*viz* v1.9.x)

The book covers the same recipe of introducing features in the opening chapters before focussing on individual topics such as: scientific programming, graphics, machine learning etc.

It assumes you already have some basic working knowledge of a high-level language such as MATLAB, R, Python and  Javascript.

This book covers the *(some of)* following features:

* Install Julia, setup the environment and add various packages.
* Develop *own* types to extend the built-in type system
* Explore, use and create macros.
* Integrate Julia with other languages such as C, Python, and Perl.
* Process data sets and apply statistical analysis. 
* Apply scientific programming techniques to solution of dynamic systems, optimisation etc.  
* Visualise your data in Julia with plotting packages.
* Interface with SQL and NoSQL databases.
* Examine Julia's approach to Bayesian methods and machine learning.


### Instructions and Navigations

All of the code is organised into folders.

For example from the first chapter.
```
function basel(M::Integer)
  @assert M > 0
  s = 0.0
  for i = 1:M
    s += 1.0/float(i)^2
  end
  return s
end
@show(basel(10^6));
basel(10^6) = 1.6449330668487
```

### Software List

With the following software  you can run all code files present in the book (Chapters 01-11).

| Chapter | Software Required | Operating System |
| --- | --- | --- |
| All | Julia v1.9.x | Windows, Linux, macOS |



### Julia's Resources
- Homepage: https://julialang.org
- Binaries: https://julialang.org/downloads/
- Source code: https://github.com/JuliaLang/julia
- Documentation: https://docs.julialang.org
- Packages: https://julialang.org/packages/
- Discussion forum: https://discourse.julialang.org


### Get to Know the Author

Malcolm Sherrington has been working in computing for over 40 years. He has degrees in mathematics, chemistry, and engineering and has held lectureships at two separate universities in UK. Additionally in the 1980's he established his own consultancy initially working in Aerospace establishments in both Rome and Munich and then returned to Britain to work on medical imaging at a major NHS research hospital and subsequently on as a full stack developer in the banking sector of the City of London.

He started programming scientific problems in the early 1970s using Fortran, as part of his graduate and then doctoral studies then went on to work on engineering systems written in DEC PDP/VAX assembly programming and as an early adopter of Unix moved via C/C++, through to Ada and Common Lisp.  Latterly he was involved with data processing and analytics in Perl, Python, and R, in particular working with large scale SQL and then NoSQL databases in financial establishments.

Malcolm was the founder of the London Julia User Group in 2012, as well as then being a co-organizer of the UK High Performance Computing and the Financial Engineers and Quant (London) meetup groups. He has spoken widely in both the UK and Europe on topics related to Julia and Financial Engineering.
