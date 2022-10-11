using ElasticArrays,Printf
using Plots

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

@views function update_old!((;τxx_old,τyy_old,τxy_old,τxx,τyy,τxy))
    τxx_old .= τxx
    τyy_old .= τyy
    τxy_old .= τxy
    return
end

@views function update_iteration_params!((;η_veτ,η_veτ_xy,dτ_ρx,dτ_ρy,Gdτ,Gdτ_xy,η,η_xy,G,G_xy,ητ),dt,re_mech,vpdτ,lτ,r)
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

@views function update_stresses!((;Pr,τxx,τyy,τxy,τxx_old,τyy_old,τxy_old,Vx,Vy,∇V,η_veτ,η_veτ_xy,η,η_xy,Gdτ,G,G_xy),r,dt,dx,dy)
    ∇V   .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
    Pr  .-= r.*Gdτ.*∇V
    τxx .+= (.-(τxx .- τxx_old)./(G.*dt) .- τxx./η .+ 2.0.*(diff(Vx,dims=1)./dx .- ∇V./3.0)).*η_veτ
    τyy .+= (.-(τyy .- τyy_old)./(G.*dt) .- τyy./η .+ 2.0.*(diff(Vy,dims=2)./dy .- ∇V./3.0)).*η_veτ
    τxy[2:end-1,2:end-1] .+= (.-(τxy[2:end-1,2:end-1] .- τxy_old[2:end-1,2:end-1])./(G_xy.*dt) .- τxy[2:end-1,2:end-1]./η_xy .+ (diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)).*η_veτ_xy
    return
end

@views function update_velocities!((;Vx,Vy,Pr,τxx,τyy,τxy,dτ_ρx,dτ_ρy,η,ρgx,ρgy),dx,dy)
    Vx[2:end-1,:] .+= dτ_ρx.*(diff(.-Pr.+τxx,dims=1)./dx .+ diff(τxy[2:end-1,:],dims=2)./dy .- ρgx)
    Vy[:,2:end-1] .+= dτ_ρy.*(diff(.-Pr.+τyy,dims=2)./dy .+ diff(τxy[:,2:end-1],dims=1)./dx .- ρgy)
    return
end

@views function compute_residuals!((;r_Vx,r_Vy,Pr,τxx,τyy,τxy,ρgx,ρgy),dx,dy)
    r_Vx .= diff(.-Pr[:,2:end-1].+τxx[:,2:end-1],dims=1)./dx .+ diff(τxy[2:end-1,2:end-1],dims=2)./dy .- ρgx[:,2:end-1]
    r_Vy .= diff(.-Pr[2:end-1,:].+τyy[2:end-1,:],dims=2)./dy .+ diff(τxy[2:end-1,2:end-1],dims=1)./dx .- ρgy[2:end-1,:]
    return
end

function main()
    # physics
    lx, ly      = 10.0, 10.0
    η0          = 1.0
    ηi          = 0.1
    G0          = 1.0
    ρg0         = 1.0
    ξ           = 1.0         # Maxwell relaxation time
    rad         = 0.1lx
    x0, y0      = 0.0, 0.5*ly
    dt          = 1.0
    # numerics
    nx, ny      = 127, 127
    nt          = 1
    ϵtol        = (1e-6, 1e-6, 1e-6)
    maxiter     = 100max(nx,ny)
    ncheck      = ceil(Int,5max(nx,ny))
    r           = 0.7
    re_mech     = 2π
    # preprocessing
    dt          = η0/(G0*ξ + 1e-15)
    dx,dy       = lx/nx,ly/ny
    xv,yv       = LinRange(-lx/2,lx/2,nx+1),LinRange(0,ly,ny+1)
    xc,yc       = amean1(xv),amean1(yv)
    lτ          = min(lx,ly)
    vpdτ        = min(dx,dy)/sqrt(2.1)
    # array allocation
    fields = (
        Vx       = zeros(nx+1,ny  ),
        Vy       = zeros(nx  ,ny+1),
        Pr       = zeros(nx  ,ny  ),
        ∇V       = zeros(nx  ,ny  ),
        τxx      = zeros(nx  ,ny  ),
        τyy      = zeros(nx  ,ny  ),
        τxy      = zeros(nx+1,ny+1),
        τxx_old  = zeros(nx  ,ny  ),
        τyy_old  = zeros(nx  ,ny  ),
        τxy_old  = zeros(nx+1,ny+1),
        η        = zeros(nx  ,ny  ),
        η_xy     = zeros(nx-1,ny-1),
        τII      = zeros(nx  ,ny  ),
        Vmag     = zeros(nx  ,ny  ),
        η_veτ    = zeros(nx  ,ny  ),
        η_veτ_xy = zeros(nx-1,ny-1),
        dτ_ρx    = zeros(nx-1,ny  ),
        dτ_ρy    = zeros(nx  ,ny-1),
        Gdτ      = zeros(nx  ,ny  ),
        Gdτ_xy   = zeros(nx-1,ny-1),
        G        = zeros(nx  ,ny  ),
        G_xy     = zeros(nx-1,ny-1),
        r_Vx     = zeros(nx-1,ny-2),
        r_Vy     = zeros(nx-2,ny-1),
        ρgy_c    = zeros(nx  ,ny  ),
        ρgx      = zeros(nx-1,ny  ),
        ρgy      = zeros(nx  ,ny-1),
        ητ       = zeros(nx  ,ny  ),
    )
    # initialisation
    fields.η .= η0
    fields.G .= G0
    rad = (xc.-x0).^2 .+ (yc'.-y0).^2
    fields.η[rad.<1] .= ηi
    fields.ρgy_c[rad.<1] .= ρg0
    fields.ρgy  .= ameany(fields.ρgy_c)
    fields.η_xy .= hmeanxy(fields.η)
    fields.G_xy .= hmeanxy(fields.G)
    iter_evo = Float64[]
    errs_evo = ElasticMatrix{Float64}(undef,length(ϵtol),0)
    opts = (aspect_ratio=1, xlims=extrema(xc), ylims=extrema(yc), c=:turbo)
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
                compute_residuals!(fields,dx,dy)
                errs = ( maximum(abs.(fields.r_Vx)),maximum(abs.(fields.r_Vy)),maximum(abs.(fields.∇V)) )
                push!(iter_evo,iter/max(nx,ny));append!(errs_evo,errs)
                @printf("  iter/nx=%1.1f,errs=[ %1.3e, %1.3e, %1.3e ]\n",iter/max(nx,ny),errs...)
            end
            iter += 1
        end
        # visualise
        fields.Vmag .= sqrt.(ameanx(fields.Vx).^2 + ameany(fields.Vy).^2)
        fields.τII  .= sqrt.(0.5.*(fields.τxx.^2 .+ fields.τyy.^2) .+ ameanxy(fields.τxy).^2)
        p1=heatmap(xc,yc,log10.(fields.η)',title="log10(η)"; opts...)
        p2=heatmap(xc,yc,fields.Pr',title="Pressure"; opts...)
        p3=heatmap(xc,yc,fields.Vmag',title="Velocity"; opts...)
        p4=heatmap(xc,yc,fields.τII',title="τII"; opts...)
        display(plot(p1,p2,p3,p4,layout=(2,2)))
    end
    return
end

main()
