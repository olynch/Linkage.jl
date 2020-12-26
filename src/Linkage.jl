module Linkage

export Lkg, lkg_model, tweak_pos, draw_linkage, interact_linkage

using Catlab.Present
using Catlab.CategoricalAlgebra.CSets
using Catlab.Graphs.BasicGraphs: TheoryGraph

using Gtk, Graphics

using JuMP, Ipopt

@present TheoryLinkage <: TheoryGraph begin
  Length :: Data
  Fixed :: Data
  length :: Attr(E,Length)
  fixed :: Attr(V,Fixed)
end

const Lkg = ACSetType(TheoryLinkage, index=[:src,:tgt]){Float64, Bool}

function lkg_model(l::Lkg,init_ps::Array{T,2}) where {T}
  model = Model()
  n = nparts(l,:V)
  @variable(model, ps[1:n,1:2])
  @variable(model, p[1:2])
  for e in 1:nparts(l,:E)
    s,t,d = subpart(l,e,:src), subpart(l,e,:tgt), subpart(l,e,:length)
    if !subpart(l,s,:fixed) || !subpart(l,t,:fixed)
      @constraint(model, sum((ps[s,:] .- ps[t,:]) .^ 2) == d^2)
    end
  end
  for v in 1:nparts(l,:V)
    if subpart(l,v,:fixed)
      @constraint(model, ps[v,:] .== init_ps[v,:])
    end
  end
  @constraint(model,p_con,p[:] .== [0.,0])
  set_optimizer(model, Ipopt.Optimizer)
  set_optimizer_attribute(model, "print_level", 0)
  model
end

function tweak_pos(model,input_ps,input_p,i)
  ps = model[:ps]
  p = model[:p]
  for i in eachindex(input_p)
    set_normalized_rhs(model[:p_con][i], input_p[i])
    set_start_value(p[i], input_p[i])
  end
  for i in eachindex(input_ps)
    set_start_value(ps[i], input_ps[i])
  end
  @objective(model, Min, sum((ps[i,:] .- p[:]) .^ 2))
  optimize!(model)
  value.(ps)
end

radius = 20

function draw_linkage(ctx,l,ps)
  set_source_rgb(ctx, 0, 0, 0)
  for v in 1:nparts(l,:V)
    arc(ctx, ps[v,1], ps[v,2], radius, 0, 2pi)
    stroke(ctx)
  end
  for e in 1:nparts(l,:E)
    s,t = subpart(l,e,:src), subpart(l,e,:tgt)
    move_to(ctx,ps[s,1],ps[s,2])
    line_to(ctx,ps[t,1],ps[t,2])
    stroke(ctx)
  end
end

function draw_path(ctx, path)
  if length(path) > 0
    set_source_rgb(ctx, 1, 0, 0)
    move_to(ctx, path[1][1], path[1][2])
    for pt in path[2:end]
        line_to(ctx, pt[1], pt[2])
    end
    stroke(ctx)
  end
end

peaucellierLkg = @acset Lkg begin
  V = 6
  E = 7
  src = [2,3,3,6,5,1,1]
  tgt = [3,6,5,4,4,5,6]
  length = [100.,100,100,100,100,280,280]
  fixed = [true,true,false,false,false,false]
end

peaucellierPs = [300. 300 ; 400 300 ; 500 300 ; 642 300 ; 571 250 ; 571 350 ]

wattsLkg = @acset Lkg begin
  V = 5
  E = 5
  src = [1,2,3,2,3]
  tgt = [2,3,4,5,5]
  length = [500,100,500,50,50]
  fixed = [true,false,false,true,false]
end

wattsPs = [200. 300 ; 700 350 ; 700 250 ; 1200 300 ; 700 300 ]

function interact_linkage(l,init_ps)
  n = nparts(l,:V)
  c = @GtkCanvas()
  ps = copy(init_ps)
  model = lkg_model(l,init_ps)
  selected = nothing
  win = GtkWindow(c,"Linkage Interaction")
  cur_path = Vector{Float64}[]
  @guarded draw(c) do widget
    ctx = getgc(c)
    w = width(c)
    h = height(c)
    set_source_rgb(ctx,1,1,1)
    rectangle(ctx,0,0,w,h)
    fill(ctx)
    draw_path(ctx, cur_path)
    draw_linkage(ctx,l,ps[:,:])
  end
  c.mouse.button1press = @guarded (widget, event) -> begin
    selected = findfirst(i -> sum((ps[i,:] .- [event.x, event.y]) .^ 2) < radius^2, 1:n)
    cur_path = [ps[selected,:]]
  end
  c.mouse.button1release = @guarded (widget, event) -> begin
    selected = nothing
    cur_path = []
    draw(c,true)
  end
  c.mouse.motion = @guarded (widget, event) -> begin
    if selected != nothing
      ps[:,:] = tweak_pos(model,ps,Float64[event.x,event.y],selected)
      push!(cur_path, ps[selected,:])
      draw(c, true)
    end
  end
  show(c)
end

end
