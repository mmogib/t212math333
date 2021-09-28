### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ be68235e-139c-11ec-150d-69521e016a4a
using PlutoUI

# ╔═╡ 096ad27c-834a-49fb-a0c1-bcb8beaefffb
begin
	y(n::Int64) = n == 1 ? 1 : (1/4)*(2*y(n-1)+3)
	y(1)
end

# ╔═╡ 0eadad40-6dfb-4ca9-8af4-a6c5fa77f341
Slder = @bind slder Slider(1:10, show_value=true) 

# ╔═╡ 4a950791-b6eb-466d-8fda-35de756c469b
begin
	s(n::Int64,a::T) where {T<:Real} = n == 1 ? 1 : 0.5*(s(n-1,a)+a/s(n-1,a))
	s(slder,2)
end

# ╔═╡ Cell order:
# ╠═096ad27c-834a-49fb-a0c1-bcb8beaefffb
# ╠═0eadad40-6dfb-4ca9-8af4-a6c5fa77f341
# ╠═4a950791-b6eb-466d-8fda-35de756c469b
# ╠═be68235e-139c-11ec-150d-69521e016a4a
