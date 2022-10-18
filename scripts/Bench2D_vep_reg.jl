using ElasticArrays,Printf
using Plots,Plots.Measures
default(size=(800,500),framestyle=:box,label=false,grid=false,margin=3mm,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

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
@views maxloc(A)  = max.(A[1:end-2,1:end-2],A[1:end-2,2:end-1],A[1:end-2,3:end],
                         A[2:end-1,1:end-2],A[2:end-1,2:end-1],A[2:end-1,3:end],
                         A[3:end  ,1:end-2],A[3:end  ,2:end-1],A[3:end  ,3:end])
@views bc2!(A)    = begin A[[1,end],:]=A[[2,end-1],:]; A[:,[1,end]]=A[:,[2,end-1]]; end
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
    ητ[2:end-1,2:end-1] .= maxloc(η)
    bc2!(ητ)
    Gdτ      .= ητ.*(re_mech/lτ*vpdτ/(r+2.0))
    Gdτ_xy   .= avxy(Gdτ)
    dτ_ρx    .= vpdτ*lτ./re_mech./avx(ητ)
    dτ_ρy    .= vpdτ*lτ./re_mech./avy(ητ)
    η_veτ    .= 1.0./(1.0./Gdτ    .+ 1.0./η    .+ 1.0./(G.*dt))
    η_veτ_xy .= 1.0./(1.0./Gdτ_xy .+ 1.0./η_xy .+ 1.0./(G_xy.*dt))
    return
end

@views function update_stresses!((;εxx1,εyy1,εxy1,εII1,Pr,εxx,εyy,εxy,εxyv,εII,dτxx,dτyy,dτxy,τxx,τyy,τxy,τxyv,τxx_old,τyy_old,τxy_old,Vx,Vy,∇V,η_veτ,η_veτ_xy,η,η_xy,Gdτ,G,G_xy,ηb,F,λ,dQdτxx,dQdτyy,dQdτxy,τII,η_vep,Fchk),τ_y,sinϕ,η_reg,r,dt,dx,dy)
    ∇V   .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
    Pr  .-= r.*Gdτ.*∇V
    # τxx .+= (.-(τxx .- τxx_old)./(G.*dt) .- τxx./η .+ 2.0.*(diff(Vx,dims=1)./dx .- ∇V./3.0)).*η_veτ
    # τyy .+= (.-(τyy .- τyy_old)./(G.*dt) .- τyy./η .+ 2.0.*(diff(Vy,dims=2)./dy .- ∇V./3.0)).*η_veτ
    # τxy[2:end-1,2:end-1] .+= (.-(τxy[2:end-1,2:end-1] .- τxy_old[2:end-1,2:end-1])./(G_xy.*dt) .- τxy[2:end-1,2:end-1]./η_xy .+ (diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)).*η_veτ_xy
    # strain rates
    εxx    .= diff(Vx,dims=1)./dx .- ∇V./3.0
    εyy    .= diff(Vy,dims=2)./dy .- ∇V./3.0
    εxyv[2:end-1,2:end-1] .= 0.5*(diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)
    εxy    .= ameanxy(εxyv)
    εII    .= sqrt.(0.5.*(εxx.^2 .+ εyy.^2) .+ εxy.^2)

    εxx1   .= εxx .+ 0.5.*τxx_old./(G.*dt)
    εyy1   .= εyy .+ 0.5.*τyy_old./(G.*dt)
    εxy1   .= εxy .+ 0.5.*τxy_old./(G.*dt)
    εII1   .= sqrt.(0.5.*(εxx1.^2 .+ εyy1.^2) .+ εxy1.^2)

    # stress increments
    dτxx   .= (.-(τxx .- τxx_old)./(G.*dt) .- τxx./η .+ 2.0.*εxx).*η_veτ
    dτyy   .= (.-(τyy .- τyy_old)./(G.*dt) .- τyy./η .+ 2.0.*εyy).*η_veτ
    dτxy   .= (.-(τxy .- τxy_old)./(G.*dt) .- τxy./η .+ 2.0.*εxy).*η_veτ
    τII    .= sqrt.(0.5.*((τxx.+dτxx).^2 .+ (τyy.+dτyy).^2) .+ (τxy.+dτxy).^2)
    # yield function
    F      .= τII .- τ_y .- Pr.*sinϕ
    λ      .= max.(F,0.0)./(η_veτ .+ η_reg)
    dQdτxx .= 0.5.*(τxx.+dτxx)./τII
    dQdτyy .= 0.5.*(τyy.+dτyy)./τII
    dQdτxy .=      (τxy.+dτxy)./τII
    τxx   .+= (.-(τxx .- τxx_old)./(G.*dt) .- τxx./η .+ 2.0.*(εxx .-      λ.*dQdτxx)).*η_veτ
    τyy   .+= (.-(τyy .- τyy_old)./(G.*dt) .- τyy./η .+ 2.0.*(εyy .-      λ.*dQdτyy)).*η_veτ
    τxy   .+= (.-(τxy .- τxy_old)./(G.*dt) .- τxy./η .+ 2.0.*(εxy .- 0.5.*λ.*dQdτxy)).*η_veτ
    τxyv[2:end-1,2:end-1] .= ameanxy(τxy)
    τII   .= sqrt.(0.5.*(τxx.^2 .+ τyy.^2) .+ τxy.^2)
    Fchk  .= τII .- τ_y .- Pr.*sinϕ .- λ.*η_reg
    η_vep .= τII ./ 2.0 ./ εII1
    return
end

@views function update_velocities!((;Vx,Vy,Pr,τxx,τyy,τxyv,dτ_ρx,dτ_ρy,η),dx,dy)
    Vx[2:end-1,:] .+= dτ_ρx.*(diff(.-Pr.+τxx,dims=1)./dx .+ diff(τxyv[2:end-1,:],dims=2)./dy)
    Vy[:,2:end-1] .+= dτ_ρy.*(diff(.-Pr.+τyy,dims=2)./dy .+ diff(τxyv[:,2:end-1],dims=1)./dx)
    return
end

@views function compute_residuals!((;r_Vx,r_Vy,Pr,τxx,τyy,τxyv),dx,dy)
    r_Vx .= diff(.-Pr[:,2:end-1].+τxx[:,2:end-1],dims=1)./dx .+ diff(τxyv[2:end-1,2:end-1],dims=2)./dy
    r_Vy .= diff(.-Pr[2:end-1,:].+τyy[2:end-1,:],dims=2)./dy .+ diff(τxyv[2:end-1,2:end-1],dims=1)./dx
    return
end

function main()
    # physics
    lx,ly      = 1.0,1.0
    radi       = 0.01*lx
    τ_y        = 1.6
    sinϕ       = sind(30)
    η0         = 1.0
    G0         = 1.0
    Gi         = G0/2
    ξ          = 4.0
    εbg        = 1.0
    dt         = η0/G0/ξ
    # numerics
    nx,ny      = 63,63
    nt         = 15
    η_reg      = 8.0e-3             # regularisation "viscosity"
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
    vpdτ       = 0.99*min(dx,dy)/sqrt(2.1)
    # array allocation
    fields = (
    Vx         = zeros(nx+1,ny  ),
    Vy         = zeros(nx  ,ny+1),
    Pr         = zeros(nx  ,ny  ),
    ∇V         = zeros(nx  ,ny  ),
    τxx        = zeros(nx  ,ny  ),
    τyy        = zeros(nx  ,ny  ),
    τxy        = zeros(nx  ,ny  ),
    τxyv       = zeros(nx+1,ny+1),
    τxx_old    = zeros(nx  ,ny  ),
    τyy_old    = zeros(nx  ,ny  ),
    τxy_old    = zeros(nx  ,ny  ),
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
    F          = zeros(nx  ,ny  ),
    λ          = zeros(nx  ,ny  ),
    dQdτxx     = zeros(nx  ,ny  ),
    dQdτyy     = zeros(nx  ,ny  ),
    dQdτxy     = zeros(nx  ,ny  ),
    Fchk       = zeros(nx  ,ny  ),
    εxx        = zeros(nx  ,ny  ),
    εyy        = zeros(nx  ,ny  ),
    εxy        = zeros(nx  ,ny  ),
    εII        = zeros(nx  ,ny  ),
    εxx1       = zeros(nx  ,ny  ),
    εyy1       = zeros(nx  ,ny  ),
    εxy1       = zeros(nx  ,ny  ),
    εII1       = zeros(nx  ,ny  ),
    εxyv       = zeros(nx+1,ny+1),
    dτxx       = zeros(nx  ,ny  ),
    dτyy       = zeros(nx  ,ny  ),
    dτxy       = zeros(nx  ,ny  ),
    η          = η0.*ones(nx  ,ny  ),
    η_xy       = η0.*ones(nx-1,ny-1),
    G          = G0.*ones(nx  ,ny  ),
    G_xy       = G0.*ones(nx-1,ny-1),
    η_vep      = η0.*ones(nx  ,ny  ),
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
            update_stresses!(fields,τ_y,sinϕ,η_reg,r,dt,dx,dy)
            update_velocities!(fields,dx,dy)
            if iter % ncheck == 0
                # update residuals
                compute_residuals!(fields,dx,dy)
                errs = maximum.((abs.(fields.r_Vx),abs.(fields.r_Vy),abs.(fields.∇V)))
                push!(iter_evo,iter/max(nx,ny));append!(errs_evo,errs)
                @printf("  iter/nx=%.3f,errs=[ %1.3e, %1.3e, %1.3e ] (Fchk=%1.2e)\n",iter/max(nx,ny),errs...,maximum(fields.Fchk))
            end
            iter += 1
        end
        t += dt
        push!(evo_t,t); push!(evo_τxx,maximum(fields.τxx))
        # visualisation
        fields.Vmag .= sqrt.(ameanx(fields.Vx).^2 + ameany(fields.Vy).^2)
        p1=heatmap(xc,yc,ameanx(fields.Vx)',title="Vx";opts...)
        p2=heatmap(xc,yc,fields.η_vep',title="η_vep";opts...)
        p3=heatmap(xc,yc,fields.τII',title="τII";opts...)
        p4=plot(evo_t,evo_τxx,legend=false,xlabel="time",ylabel="max(τxx)",linewidth=0,markershape=:circle,markersize=3,framestyle=:box)
        display(plot(p1,p2,p3,p4,layout=(2,2)))
    end
    return
end

main()
