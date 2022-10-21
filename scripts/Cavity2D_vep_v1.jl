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
@views maxloc(A)  = max.(A[1:end-2,1:end-2],A[1:end-2,2:end-1],A[1:end-2,3:end],
                         A[2:end-1,1:end-2],A[2:end-1,2:end-1],A[2:end-1,3:end],
                         A[3:end  ,1:end-2],A[3:end  ,2:end-1],A[3:end  ,3:end])
@views bc2!(A)    = begin A[[1,end],:]=A[[2,end-1],:]; A[:,[1,end]]=A[:,[2,end-1]]; end

@views function update_old!((;τxx_old,τyy_old,τxy_old,τxx,τyy,τxy))
    τxx_old .= τxx
    τyy_old .= τyy
    τxy_old .= τxy
    return
end

@views function update_iteration_params!((;η,ητ))
    ητ[2:end-1,2:end-1] .= maxloc(η); bc2!(ητ)
    return
end

@views function update_stresses!((;Pr,dPr,τxx,τyy,τxy,τxyv,τxx_old,τyy_old,τxy_old,εxx,εyy,εxy,εxyv,Vx,Vy,∇V,η,G,dτ_r),dt,re_mech,vdτ,lτ,r,dx,dy)
    θ_dτ   = lτ*(r+2.0)/(re_mech*vdτ)
    dτ_r  .= 1.0./(θ_dτ .+ η./(G.*dt) .+ 1.0)
    ∇V    .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
    dPr   .= .-∇V
    Pr   .+= (r/θ_dτ).*η.*dPr
    εxx   .= diff(Vx,dims=1)./dx .- ∇V./3.0
    εyy   .= diff(Vy,dims=2)./dy .- ∇V./3.0
    εxyv[2:end-1,2:end-1] .= 0.5*(diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)
    εxy   .= ameanxy(εxyv)
    τxx  .+= (.-(τxx .- τxx_old).*η./(G.*dt) .- τxx .+ 2.0.*η.*εxx).*dτ_r
    τyy  .+= (.-(τyy .- τyy_old).*η./(G.*dt) .- τyy .+ 2.0.*η.*εyy).*dτ_r
    τxy  .+= (.-(τxy .- τxy_old).*η./(G.*dt) .- τxy .+ 2.0.*η.*εxy).*dτ_r
    τxyv[2:end-1,2:end-1] .= ameanxy(τxy)
    return
end

@views function update_velocities!((;Vx,Vy,Pr,τxx,τyy,τxyv,η,ητ,ρgx,ρgy,phase),vdτ,lτ,re_mech,dx,dy)
    nudτ = vdτ*lτ/re_mech
    Vx[2:end-1,:] .+= (diff(.-Pr.+τxx,dims=1)./dx .+ diff(τxyv[2:end-1,:],dims=2)./dy .- ρgx).*nudτ./avx(ητ)
    Vy[:,2:end-1] .+= (diff(.-Pr.+τyy,dims=2)./dy .+ diff(τxyv[:,2:end-1],dims=1)./dx .- ρgy).*nudτ./avy(ητ)
    Vx[:,1]       .= 1.0./3.0.*Vx[:,2]
    # Vx[end,:]     .= Vx[end-1,:] .+ 0.5.*dx./η[end,:].*(Pr[end,:] .+ 1.0./3.0.*(-Pr[end-1,:] .+ 2.0.*η[end-1,:].*(Vx[end-1,:] .- Vx[end-2,:])./dx))
    Vy[:,end]     .= Vy[:,end-1] .+ 0.5.*dy./η[:,end].*(Pr[:,end] .+ 1.0./3.0.*(-Pr[:,end-1] .+ 2.0.*η[:,end-1].*(Vy[:,end-1] .- Vy[:,end-2])./dy))
    Vy[:,1]       .= (1 .- phase[:,1]) .* (Vy[:,2] .- dy ./ (2.0 .* η[:,1]) .* (Pr[:,1] .+ 1.0./3.0.*(-Pr[:,2] .+ 2.0.*η[:,2].*(Vy[:,3] .- Vy[:,2])./dy)))
    return
end

@views function compute_residuals!((;r_Vx,r_Vy,Pr,τxx,τyy,τxyv,ρgx,ρgy),dx,dy)
    r_Vx .= diff(.-Pr[:,2:end-1].+τxx[:,2:end-1],dims=1)./dx .+ diff(τxyv[2:end-1,2:end-1],dims=2)./dy .- ρgx[:,2:end-1]
    r_Vy .= diff(.-Pr[2:end-1,:].+τyy[2:end-1,:],dims=2)./dy .+ diff(τxyv[2:end-1,2:end-1],dims=1)./dx .- ρgy[2:end-1,:]
    return
end

@views function compte_η_G_ρg!((;η,G,ρgy_c,phase,ηb),η0,G0,ρg0,xc,yc,x0,y0c,y0d,r_cav,r_dep)
    Threads.@threads for iy in axes(η,2)
        for ix in axes(η,1)
            sd_air = min(sqrt((xc[ix]-x0)^2 + 2*(yc[iy]-y0c)^2)-r_cav,
                         sqrt((xc[ix]-x0)^2 + 5*(yc[iy]-y0d)^2)-r_dep)
            t_air  = 0.5*(tanh(-sd_air/0.06) + 1)
            t_ice  = 1.0 - t_air
            η[ix,iy]     = t_ice*η0.ice  + t_air*η0.air
            G[ix,iy]     = t_ice*G0.ice  + t_air*G0.air
            ρgy_c[ix,iy] = t_ice*ρg0.ice + t_air*ρg0.air
            phase[ix,iy] = 1.0 - t_air
            ηb[ix,iy]    = (1.0 - t_air)*1e12 + t_air*1.0
        end
    end
    return
end

function main()
    # physics
    lx,ly      = 20.0,10.0
    η0         = (ice = 1.0 , air = 1e-6)
    G0         = (ice = 1.0 , air = 1e6 )
    ρg0        = (ice = 0.9 , air = 0.0 )
    r_cav      = 0.4*min(lx,ly)
    r_dep      = 1.5*min(lx,ly)
    x0,y0c,y0d = 0.0,0.0*ly,1.4*ly
    ξ          = 1.0
    dt         = η0.ice/(G0.ice*ξ + 1e-15)
    # numerics
    nx         = 256
    ny         = ceil(Int,nx*ly/lx)
    nt         = 1
    ϵtol       = (1e-6,1e-6,1e-6)
    maxiter    = 50max(nx,ny)
    ncheck     = ceil(Int,5max(nx,ny))
    r          = 0.7
    re_mech    = 3π
    # preprocessing
    dx,dy      = lx/nx,ly/ny
    xv,yv      = LinRange(-lx/2,lx/2,nx+1),LinRange(0,ly,ny+1)
    xc,yc      = amean1(xv),amean1(yv)
    lτ         = min(lx,ly)
    vdτ        = 0.5*min(dx,dy)/sqrt(2.1)
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
    εxx        = zeros(nx  ,ny  ),
    εyy        = zeros(nx  ,ny  ),
    εxy        = zeros(nx  ,ny  ),
    εxyv       = zeros(nx+1,ny+1),
    η          = zeros(nx  ,ny  ),
    τII        = zeros(nx  ,ny  ),
    Vmag       = zeros(nx  ,ny  ),
    dPr        = zeros(nx  ,ny  ),
    G          = zeros(nx  ,ny  ),
    r_Vx       = zeros(nx-1,ny-2),
    r_Vy       = zeros(nx-2,ny-1),
    dτ_r       = zeros(nx  ,ny  ),
    ρgy_c      = zeros(nx  ,ny  ),
    ρgx        = zeros(nx-1,ny  ),
    ρgy        = zeros(nx  ,ny-1),
    phase      = zeros(nx  ,ny  ),
    ηb         = zeros(nx  ,ny  ),
    ητ         = zeros(nx  ,ny  ),
    )
    # initialisation
    compte_η_G_ρg!(fields,η0,G0,ρg0,xc,yc,x0,y0c,y0d,r_cav,r_dep)
    fields.Pr   .= reverse(cumsum(reverse(fields.ρgy_c,dims=2),dims=2).*dy,dims=2)
    fields.ρgy  .= ameany(fields.ρgy_c)
    iter_evo = Float64[]
    errs_evo = ElasticMatrix{Float64}(undef,length(ϵtol),0)
    opts = (aspect_ratio=1, xlims=extrema(xc), ylims=extrema(yc), c=:turbo, framestyle=:box)
    mask = copy(fields.phase); mask[mask.<0.7].=NaN
    t = 0.0; evo_t=[]; evo_τxx=[]
    # time loop
    for it = 1:nt
        @printf("it=%d\n",it)
        update_old!(fields)
        errs = 2.0.*ϵtol; iter = 1
        resize!(iter_evo,0); resize!(errs_evo,length(ϵtol),0)
        while any(errs .>= ϵtol) && iter <= maxiter
            update_iteration_params!(fields)
            update_stresses!(fields,dt,re_mech,vdτ,lτ,r,dx,dy)
            update_velocities!(fields,vdτ,lτ,re_mech,dx,dy)
            if iter % ncheck == 0
                # update residuals
                compute_residuals!(fields,dx,dy)
                errs = maximum.((abs.(fields.r_Vx),abs.(fields.r_Vy),abs.(fields.dPr)))
                push!(iter_evo,iter/max(nx,ny));append!(errs_evo,errs)
                @printf("  iter/nx=%.3f,errs=[ %1.3e, %1.3e, %1.3e ]\n",iter/max(nx,ny),errs...)
            end
            iter += 1
        end
        t += dt
        push!(evo_t,t); push!(evo_τxx,maximum(fields.τxx))
        # visualisation
        fields.Vmag .= sqrt.(ameanx(fields.Vx).^2 + ameany(fields.Vy).^2)
        fields.τII  .= sqrt.(0.5.*(fields.τxx.^2 .+ fields.τyy.^2) .+ fields.τxy.^2)
        p1=heatmap(xc,yc,log10.(fields.η)',title="log10(η)";opts...)
        # p2=heatmap(xc,yc,mask' .* fields.Pr',title="Pressure";opts...)
        p2=heatmap(xc,yc,mask' .* fields.τII',title="Pressure";opts...)
        p3=heatmap(xc,yc,mask' .* fields.Vmag',title="Vmag";opts...)
        p4=plot(evo_t,evo_τxx,legend=false,xlabel="time",ylabel="max(τxx)",linewidth=0,markershape=:circle,markersize=3,framestyle=:box)
        display(plot(p1,p2,p3,p4,layout=(2,2)))
    end
    return
end

main()
