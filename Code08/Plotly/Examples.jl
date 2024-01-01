# [ A box plot of stock prices ]
#=
  Stock market values for Google, Apple, Amazon, and Microsoft 
  for the years 2018-19. All stocks are normalized to start with 
  their closing values on 01-01-2019.
  
cd(ENV["HOME"]*"/MJ2")

using DataFrames, CSV, PlotlyJS
df = CSV.read("Datasources/stocks4.csv",DataFrame);

trace1 = box(;x = df[!,:GOOG], name="Google");
trace2 = box(;x = df[!,:AAPL], name="Apple");
trace3 = box(;x = df[!,:AMZN], name="Amazon");
trace4 = box(;x = df[!,:MSFT], name="Microsoft");

data = [trace1, trace2, trace3, trace4];
layout = Layout(;title="Selected Stock Values for 2018-2019",
                xaxis=attr(title="Normalized to 1.0 on 1st January 2018",
                zeroline=false));
plot(data, layout)

## Using Plotly [online]
#=
  The API key is (currently) a 20-character ASCII string returned on registration 
  by Plotly and can be viewed, if forgotten, by logging on to https://Plot.ly
=#
Plotly.set_credentials_file(Dict("username"=>"sherrinm",
                                 "api_key"=>"<<API key>") )
# Sign in to Plotly
using Plotly
Plotly.signin("sherrinm", "<<Plotly API key>>")

# A contour plot of some sinusoids with a randomly generated component
N = 100;
X = collect(range(-2*pi, stop=2*pi, length=N));
Y = collect(range(-2*pi, stop=2*pi, length=N));
Z = zeros(N, N);

for i = 1:N, j = 1:N
  r2 = (X[i]^2 + Y[j]^2)
  Z[i,j] = sin(X[i]) * cos(Y[j]) * sin(r2)/log(r2+1)
end

data = contour(z=Z , x=>X,  y=Y);
resp = Plotly.plot(data);
plot_url = Plotly.post(resp);
 
