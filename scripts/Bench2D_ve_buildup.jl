using ElasticArrays,Printf
using Plots,Plots.Measures
default(size=(800,500),framestyle=:box,label=false,grid=false,margin=3mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)
# using GLMakie
# include("vis_helpers.jl")

@inline amean(a,b) = 0.5*(a + b)
@inline hmean(a,b) = 2.0/(1.0/a + 1.0/b)
@inline amean4(a,b,c,d) = 0.25*(a+b+c+d)
@inline hmean4(a,b,c,d) = 4.0/(1.0/a+1.0/b+1.0/c+1.0/d)
const av  = amean
const av4 = amean4
@views amean1(A)  = 0.5.*(A[1:end-1] .+ A[2:end])
@views avx(A)     = av.(A[1:end-1,:], A[2:end,:])
@views avy(A)     = av.(A[:,1:end-1], A[:,2:end])
@views avxy(A)    = av4.(A[1:end-1,1:end-1],A[2:end,1:end-1],A[1:end-1,2:end],A[2:end,2:end])
@views ameanx(A)  = amean.(A[1:end-1,:], A[2:end,:])
@views ameany(A)  = amean.(A[:,1:end-1], A[:,2:end])
@views ameanxy(A) = amean4.(A[1:end-1,1:end-1],A[2:end,1:end-1],A[1:end-1,2:end],A[2:end,2:end])
@views hmeanx(A)  = hmean.(A[1:end-1,:], A[2:end,:])
@views hmeany(A)  = hmean.(A[:,1:end-1], A[:,2:end])
@views hmeanxy(A) = hmean4.(A[1:end-1,1:end-1],A[2:end,1:end-1],A[1:end-1,2:end],A[2:end,2:end])

# 1. Viscous
# (τ - τ_it)/Gdτ + τ/η = 2*e
# τ*(1/Gdτ + 1/η) - τ_it/Gdτ = 2*e
# τ += (-τ/η + 2*e)/(1/Gdτ + 1/η)
# 2. Maxwell viscoelastic
# (τ - τ_it)/Gdτ + τ/η + (τ - τ_old)/(G*dt) = 2*e
# τ += (-(τ - τ_old)/(G*dt) - τ/η + 2*e)/(1/Gdτ + 1/η + 1/(G*dt))
# Pressure
# 1./Kdτ.*(Pr - Pr_it) + Pr./ηb = -∇V
# Pr .+= (.-Pr./ηb .- ∇V)./(1.0/Kdτ + 1.0./ηb)

@views function update_old!((;τxx_old,τyy_old,τxy_old,τxx,τyy,τxy))
    τxx_old .= τxx
    τyy_old .= τyy
    τxy_old .= τxy
    return
end

@views function update_iteration_params!((;η_veτ,η_veτ_xy,dτ_ρx,dτ_ρy,Gdτ,Gdτ_xy,η,η_xy,G,G_xy,ηb,ητ),dt,re_mech,vpdτ,lτ,r)
    ητ[2:end-1,2:end-1]  .= max.(η[2:end-1,2:end-1],η[1:end-2,2:end-1],η[3:end  ,2:end-1],
                                 η[2:end-1,1:end-2],η[2:end-1,3:end  ],η[1:end-2,1:end-2],
                                 η[3:end  ,1:end-2],η[1:end-2,3:end  ],η[3:end  ,3:end  ],)
    ητ[[1,end],:] .= ητ[[2,end-1],:]; ητ[:,[1,end]] .= ητ[:,[2,end-1]]
    Gdτ      .= ητ.*(re_mech/lτ*vpdτ/(r+2.0))
    Gdτ_xy   .= avxy(Gdτ)
    dτ_ρx    .= vpdτ*lτ./re_mech./avx(ητ)
    dτ_ρy    .= vpdτ*lτ./re_mech./avy(ητ)
    η_veτ    .= 1.0./(1.0./Gdτ    .+ 1.0./η    .+ 1.0./(G.*dt))
    η_veτ_xy .= 1.0./(1.0./Gdτ_xy .+ 1.0./η_xy .+ 1.0./(G_xy.*dt))
    return
end

@views function update_stresses!((;Pr,τxx,τyy,τxy,τxx_old,τyy_old,τxy_old,Vx,Vy,∇V,η_veτ,η_veτ_xy,η,η_xy,Gdτ,G,G_xy,ηb),r,dt,dx,dy)
    nx,ny = size(Pr)
    ∇V   .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
    Pr  .-= r.*Gdτ.*∇V
    τxx .+= (.-(τxx .- τxx_old)./(G.*dt) .- τxx./η .+ 2.0.*(diff(Vx,dims=1)./dx .- ∇V./3.0)).*η_veτ
    τyy .+= (.-(τyy .- τyy_old)./(G.*dt) .- τyy./η .+ 2.0.*(diff(Vy,dims=2)./dy .- ∇V./3.0)).*η_veτ
    τxy[2:end-1,2:end-1] .+= (.-(τxy[2:end-1,2:end-1] .- τxy_old[2:end-1,2:end-1])./(G_xy.*dt) .- τxy[2:end-1,2:end-1]./η_xy .+ (diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)).*η_veτ_xy
    return
end

@views function update_velocities!((;Vx,Vy,Pr,τxx,τyy,τxy,dτ_ρx,dτ_ρy,η),dx,dy)
    nx,ny = size(Pr)
    Vx[2:end-1,:] .+= dτ_ρx.*(diff(.-Pr.+τxx,dims=1)./dx .+ diff(τxy[2:end-1,:],dims=2)./dy)
    Vy[:,2:end-1] .+= dτ_ρy.*(diff(.-Pr.+τyy,dims=2)./dy .+ diff(τxy[:,2:end-1],dims=1)./dx)
    return
end

@views function compute_residuals!((;r_Vx,r_Vy,Pr,τxx,τyy,τxy),dx,dy)
    r_Vx .= diff(.-Pr[:,2:end-1].+τxx[:,2:end-1],dims=1)./dx .+ diff(τxy[2:end-1,2:end-1],dims=2)./dy
    r_Vy .= diff(.-Pr[2:end-1,:].+τyy[2:end-1,:],dims=2)./dy .+ diff(τxy[2:end-1,2:end-1],dims=1)./dx
    return
end

function main()
    # physics
    lx,ly      = 1.0,1.0
    radi       = 0.01*lx
    η0         = 1.0
    G0         = 1.0
    Gi         = 0.5
    ξ          = 2.0
    εbg        = 1.0
    dt         = η0/(G0*ξ + 1e-15)
    # numerics
    nx,ny      = 127,127
    nt         = 10
    ϵtol       = (1e-6,1e-6,1e-6)
    maxiter    = 50max(nx,ny)
    ncheck     = ceil(Int,5max(nx,ny))
    r          = 0.7
    re_mech    = 3π
    # preprocessing
    dx,dy      = lx/nx,ly/ny
    xv,yv      = LinRange(-lx/2,lx/2,nx+1),LinRange(-ly/2,ly/2,ny+1)
    xc,yc      = amean1(xv),amean1(yv)
    lτ         = min(lx,ly)
    vpdτ       = 1.0*min(dx,dy)/sqrt(2.1)
    # array allocation
    fields = (
    Vx         = zeros(nx+1,ny  ),
    Vy         = zeros(nx  ,ny+1),
    Pr         = zeros(nx  ,ny  ),
    ∇V         = zeros(nx  ,ny  ),
    τxx        = zeros(nx  ,ny  ),
    τyy        = zeros(nx  ,ny  ),
    τxy        = zeros(nx+1,ny+1),
    τxx_old    = zeros(nx  ,ny  ),
    τyy_old    = zeros(nx  ,ny  ),
    τxy_old    = zeros(nx+1,ny+1),
    τII        = zeros(nx  ,ny  ),
    Vmag       = zeros(nx  ,ny  ),
    η_veτ      = zeros(nx  ,ny  ),
    η_veτ_xy   = zeros(nx-1,ny-1),
    dτ_ρx      = zeros(nx-1,ny  ),
    dτ_ρy      = zeros(nx  ,ny-1),
    Gdτ        = zeros(nx  ,ny  ),
    Gdτ_xy     = zeros(nx-1,ny-1),
    r_Vx       = zeros(nx-1,ny-2),
    r_Vy       = zeros(nx-2,ny-1),
    ηb         = zeros(nx  ,ny  ),
    ητ         = zeros(nx  ,ny  ),
    η          = η0.*ones(nx  ,ny  ),
    η_xy       = η0.*ones(nx-1,ny-1),
    G          = G0.*ones(nx  ,ny  ),
    G_xy       = G0.*ones(nx-1,ny-1),
    )
    # initialisation
    (Xvx,Yvx) = ([x for x=xv,y=yc], [y for x=xv,y=yc])
    (Xvy,Yvy) = ([x for x=xc,y=yv], [y for x=xc,y=yv])
    rad       = xc.^2 .+ yc'.^2
    fields.G[rad.<radi] .= Gi
    fields.Vx   .=   εbg.*Xvx
    fields.Vy   .= .-εbg.*Yvy
    fields.η_xy .= hmeanxy(fields.η)
    fields.G_xy .= hmeanxy(fields.G)
    iter_evo = Float64[]; errs_evo = ElasticMatrix{Float64}(undef,length(ϵtol),0)
    opts = (aspect_ratio=1, xlims=extrema(xc), ylims=extrema(yc), c=:turbo, framestyle=:box)
    t = 0.0; evo_t=[]; evo_τxx=[]
    # time loop
    for it = 1:nt
        @printf("it=%d\n",it)
        update_old!(fields)
        errs = 2.0.*ϵtol; iter = 1
        resize!(iter_evo,0); resize!(errs_evo,length(ϵtol),0)
        while any(errs .>= ϵtol) && iter <= maxiter
            update_iteration_params!(fields,dt,re_mech,vpdτ,lτ,r)
            update_stresses!(fields,r,dt,dx,dy)
            update_velocities!(fields,dx,dy)
            if iter % ncheck == 0
                # update residuals
                compute_residuals!(fields,dx,dy)
                errs = maximum.((abs.(fields.r_Vx),abs.(fields.r_Vy),abs.(fields.∇V)))
                push!(iter_evo,iter/max(nx,ny));append!(errs_evo,errs)
                @printf("  iter/nx=%.3f,errs=[ %1.3e, %1.3e, %1.3e ]\n",iter/max(nx,ny),errs...)
            end
            iter += 1
        end
        t += dt
        push!(evo_t,t); push!(evo_τxx,maximum(fields.τxx))
        # visualisation
        fields.Vmag .= sqrt.(ameanx(fields.Vx).^2 + ameany(fields.Vy).^2)
        fields.τII  .= sqrt.(0.5.*(fields.τxx.^2 .+ fields.τyy.^2) .+ ameanxy(fields.τxy).^2)
        p1=heatmap(xc,yc,ameanx(fields.Vx)',title="Vx";opts...)
        p2=heatmap(xc,yc,fields.τII',title="τII";opts...)
        p3=heatmap(xc,yc,fields.G',title="G";opts...)
        p4=plot(evo_t,evo_τxx,legend=false,xlabel="time",ylabel="max(τxx)",linewidth=0,markershape=:circle,markersize=3,framestyle=:box)
        display(plot(p1,p2,p3,p4,layout=(2,2)))
    end
    return
end

main()
