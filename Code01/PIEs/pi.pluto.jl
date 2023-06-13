### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 1387ac0b-dc6a-4d1a-ac5e-b1c81136efc0
sumsq(x,y) = x*x + y*y;

# ╔═╡ 53948eff-8c5c-4527-8379-e2a50f5c4672
N = 1000000

# ╔═╡ 6920e3e5-d92b-4043-918d-efaf142e9a86
begin
K = 0
for i in 1:N
 if sumsq(rand(), rand()) < 1.0
    K += 1
 end
end
end

# ╔═╡ bbac9f27-0aff-4555-ad7c-ee15291866b1
P = 4.0*(K / N)

# ╔═╡ 6f1cf3e8-9590-434c-a0c0-cbba27c832f6
md"Estimate of PI for $N trials is $P"

# ╔═╡ Cell order:
# ╠═1387ac0b-dc6a-4d1a-ac5e-b1c81136efc0
# ╠═53948eff-8c5c-4527-8379-e2a50f5c4672
# ╠═6920e3e5-d92b-4043-918d-efaf142e9a86
# ╠═bbac9f27-0aff-4555-ad7c-ee15291866b1
# ╠═6f1cf3e8-9590-434c-a0c0-cbba27c832f6
