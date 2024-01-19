# Mastering Julia

<a href="https://www.packtpub.com/product/mastering-julia-second-edition/9781805129790"><img src="https://content.packt.com/B21060/cover_image_small.jpg" alt="no-image" height="256px" align="right"></a>

This is the code repository for [Mastering Julia](https://www.packtpub.com/product/mastering-julia-second-edition/9781805129790), published by Packt.

**Enhance your analytical and programming skills for data modeling and processing with Julia**

## What is this book about?
This book covers the application of Julia v1.8.x in the areas of scientific computing and data science. The book will be of use to those with some previous knowledge of Julia or as a primer for programmers familiar with other scripting or compiled languages.

This book covers the following exciting features:
* Develop simple scripts in Julia using the REPL, code editors, and web-based IDEs
* Get to grips Julia’s type system, multiple dispatch, metaprogramming, and macro development
* Interact with data files, tables, data frames, SQL, and NoSQL databases
* Delve into statistical analytics, linear programming, and optimization problems
* Create graphics and visualizations to enhance modeling and simulation in Julia
* Understand Julia's main approaches to machine learning, Bayesian analysis, and AI

If you feel this book is for you, get your [copy](https://www.amazon.com/Mastering-Julia-analytical-programming-processing/dp/1805129791/ref=sr_1_1?crid=SORTAF1XG4TP&keywords=Mastering+Julia&qid=1705575656&sprefix=mastering+jul%2Caps%2C308&sr=8-1) today!


## Instructions and Navigations
All of the code is organized into folders. For example, folder Code02 contains code for Chapter 2.

The code will look like the following:
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


**Following is what you need for this book:**
This book is not an introduction to computer programming, but a practical guide for developers who want to enhance their basic knowledge of Julia, or those wishing to augment their skill set by adding Julia to their existing roster of programming languages. Familiarity with a scripting language such as Python or R, or a compiled language such as C/C++, C# or Java, is a prerequisite.

With the following software and hardware list you can run all code files present in the book (Chapter 1-11).
### Software and Hardware List
| Chapter | Software required | OS required |
| -------- | ------------------------------------ | ----------------------------------- |
| 1-11 | Julia 1.X | Windows, macOS, or Linux |


### Related products
* Hands-On Design Patterns and Best Practices with Julia [[Packt]](https://www.packtpub.com/product/hands-on-design-patterns-and-best-practices-with-julia/9781838648817) [[Amazon]](https://www.amazon.com/Hands-Design-Patterns-Julia-comprehensive/dp/183864881X/ref=sr_1_1?crid=1Z9LF0SE2EQFL&keywords=hands-on-design-patterns-and-best-practices-with-julia&qid=1705575998&sprefix=%2Caps%2C290&sr=8-1)

* Julia Programming Projects [[Packt]](https://www.packtpub.com/product/julia-programming-projects/9781788292740) [[Amazon]](https://www.amazon.com/Julia-1-0-Example-Adrian-Salceanu-ebook/dp/B076WWYYCC/ref=sr_1_1?crid=23YOU21CN49LX&keywords=Julia+Programming+Projects&qid=1705576041&sprefix=hands-on-design-patterns-and-best-practices-with-julia%2Caps%2C315&sr=8-1)

## Get to Know the Author
**Malcolm Sherrington**
has been working in computing for over 40 years. Malcolm holds degrees in Mathematics, Chemistry, and Engineering at the masters and doctorate level, and was employed as a lecturer at the University of Brighton and Imperial College London, in addition to working in the Aerospace at Rome and Munich. Subsequently, he ran his own consultancy in the Finance sector, with specific interests in Risk and High Performance Computing; augmenting his academic studies at Kings and Birkbeck Universities, London, in Financial Mathematics and Econometrics. Malcolm started programming scientific problems in Fortran and C, as part of his postgraduate research, then incorporating Ada and Lisp and recently became involved with data processing and analytics in Perl, Python, and R	

