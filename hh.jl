# hh.jl.. solve household problem for given prices wz[,z], abz[,z], taxes, etc.

function HH_root(sage, 
                 lambdain, pcin, Consin, ellin, dis_totin, yin, Ain, Savin,
                 gamin, rin, tauCin, win, tauWin, thetain, notretin, pin, taulin, ivin, abin,
                 nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi)
  
  # endogenous variables (have to be returned again)
  local lambdai  = lambdain
  local pci      = pcin
  local Consi    = Consin
  local elli     = ellin
  local dis_toti = dis_totin
  local yi       = yin
  local Ai       = Ain
  local Savi     = Savin
  
  # exogenous variables
  local gami     = gamin
  local ri       = rin
  local tauCi    = tauCin
  local wi       = win
  local tauWi    = tauWin
  local thetai   = thetain
  local notreti  = notretin
  local ppi      = pin # pi already defined
  local tauli    = taulin
  local ivi      = ivin
  local abi      = abin
  
  # EULER EQUATION: solve forward in age
  #lambdaz[sage,z]    = lambdain
  if sage < nagi
    for a in sage:(nagi-1)
      lambdai[a+1] = lambdai[a]/((1/(1+rhoi))*gami[a]*(1+ri[a]))
    end
  end
  
  # CONSUMPTION
  pci[sage:nagi]       = 1.0.+tauCi[sage:nagi]
  Consi[sage:nagi]     = (pci[sage:nagi].*lambdai[sage:nagi]).^(-sigmai)

  # HOURS SUPPLY
  elli[sage:nagi]     = ((wi[sage:nagi].*(1.0.-tauWi[sage:nagi]).*thetai[sage:nagi]./pci[sage:nagi].*(Consi[sage:nagi].^(-1/sigmai)))./parlv0i[sage:nagi]).^sigLi
  dis_toti[sage:nagi] = (sigLi/(1+sigLi)).*parlv0i[sage:nagi].*elli[sage:nagi].^((1+sigLi)/sigLi).-parlv1i[sage:nagi]
  
  # CONSUMPTION AND SAVINGS
  yi[sage:nagi]       = notreti[sage:nagi].*(wi[sage:nagi].*(1.0.-tauWi[sage:nagi]).*elli[sage:nagi].*thetai[sage:nagi]).+(1.0.-notreti[sage:nagi]).*(1.0.-tauWi[sage:nagi]).*ppi[sage:nagi].-tauli[sage:nagi]

  # ASSETS: solve forward in age
  Ai[1]         = 0
  
  if sage < nagi
    for a in sage:(nagi-1)
      Ai[a+1]   = (1+ri[a])*(Ai[a]+yi[a]+ivi[a]+abi[a]-pci[a]*Consi[a])
    end
  end
  Savi[sage:nagi]  = Ai[sage:nagi].+yi[sage:nagi].+ivi[sage:nagi].+abi[sage:nagi].-pci[sage:nagi].*Consi[sage:nagi]
  
  #return Savz[nag,z]
  return lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi
  
end

function HH(sage, maxiter, stol, atol,
           lambdain, pcin, Consin, ellin, dis_totin, yin, Ain, Savin,
           gamin, rin, tauCin, win, tauWin, thetain, notretin, pin, taulin, ivin, abin,
           nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi,
           HH_nonconvin)
  
  # endogenous variables (have to be returned again)
  local lambdai  = lambdain
  local pci      = pcin
  local Consi    = Consin
  local elli     = ellin
  local dis_toti = dis_totin
  local yi       = yin
  local Ai       = Ain
  local Savi     = Savin
  
  # exogenous variables
  local gami     = gamin
  local ri       = rin
  local tauCi    = tauCin
  local wi       = win
  local tauWi    = tauWin
  local thetai   = thetain
  local notreti  = notretin
  local ppi      = pin # pi already defined
  local tauli    = taulin
  local ivi      = ivin
  local abi      = abin
  
  local HH_nonconvi = HH_nonconvin
  
  lambdaz0 = 1.0 # initialization
  f0       = 1.0 # initialization
  
  err            = Inf
  iter           = 0
  trys           = 0
  stepsize       = 1e-6 # for numerical gradient
  
  lambdatrys     = [1.0,0.5,1.5,0.25,1.25,0.1,1.0]
  maxtrys        = length(lambdatrys)
  while_continue = true
  
  while (while_continue)
    
    while_continue = false
    lambdazsave    = lambdai[sage]
    
    while ((err > stol)||(abs(Savi[nagi]) > atol)) && (trys < maxtrys)
      
      trys += 1
      iterpertry = 0
      lambdaz1 = lambdazsave*lambdatrys[trys]
      
      breakwhile = false
      while (err > stol) && (iterpertry < maxiter) && (breakwhile == false)
        if iterpertry == 0 # Newton step for first iteration
          lambdai[sage] = lambdaz1+stepsize
          lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi = HH_root(sage, lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi, gami, ri, tauCi, wi, tauWi, thetai, notreti, ppi, tauli, ivi, abi, nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi)
          f2 = Savi[nagi]
          iter += 1
          
          if !isfinite(f2)
            breakwhile = true
            break
          end
          lambdai[sage] = lambdaz1
          lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi = HH_root(sage, lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi, gami, ri, tauCi, wi, tauWi, thetai, notreti, ppi, tauli, ivi, abi, nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi)
          f1 = Savi[nagi]
          iter += 1
          
          if !isfinite(f1)
            breakwhile = true
            break
          end
          lambdaz2 = lambdaz1 - f1*stepsize/(f2-f1)
          if (!isfinite(lambdaz2))||(lambdaz2<0)
            breakwhile = true
            break
          end
        else # Secant method
          lambdai[sage] = lambdaz1
          lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi = HH_root(sage, lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi, gami, ri, tauCi, wi, tauWi, thetai, notreti, ppi, tauli, ivi, abi, nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi)
          f1 = Savi[nagi]
          iter += 1
          
          if !isfinite(f1)
            breakwhile = true
            break
          end
          lambdaz2 = lambdaz1 - f1*(lambdaz1-lambdaz0)/(f1-f0)
          if (!isfinite(lambdaz2))||(lambdaz2<0)
            breakwhile = true
            break
          end
        end
        err = abs(lambdaz2-lambdaz1)
        lambdaz0 = lambdaz1
        lambdaz1 = lambdaz2
        f0       = f1
        iterpertry += 1
      end
    end
  end
  
  if abs(Savi[nagi]) > atol
    HH_nonconvi[nagi] = 1 # counter
  end
  
  return lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi, HH_nonconvi
  
end

function HHall(starttime, calibinit, scaleA,
           lambdain, pcin, Consin, ellin, dis_totin, yin, Ain, Savin,
           gamin, rin, tauCin, win, tauWin, thetain, notretin, pin, taulin, ivin, abin,
           nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi, ncohi, fagi, Av0i,
           HH_nonconvin)
  
    # endogenous variables (have to be returned again)
  local lambdai  = lambdain
  local pci      = pcin
  local Consi    = Consin
  local elli     = ellin
  local dis_toti = dis_totin
  local yi       = yin
  local Ai       = Ain
  local Savi     = Savin
  
  # exogenous variables
  local gami     = gamin
  local ri       = rin
  local tauCi    = tauCin
  local wi       = win
  local tauWi    = tauWin
  local thetai   = thetain
  local notreti  = notretin
  local ppi      = pin # pi already defined
  local tauli    = taulin
  local ivi      = ivin
  local abi      = abin
  
  local HH_nonconvi = HH_nonconvin
  
  Threads.@threads for z in starttime:ncohi
    if z <= nagi-fagi+starttime-1
      if calibinit
        Ai[:,z] = Av0i
      end
      Ai[nagi-(z-starttime),z] = Ai[nagi-(z-starttime),z]*scaleA
      lambdai[:,z], pci[:,z], Consi[:,z], elli[:,z], dis_toti[:,z], yi[:,z], Ai[:,z], Savi[:,z], HH_nonconvi[:,z] = HH(nagi-(z-starttime), 30, 1e-10, 0.1,lambdai[:,z], pci[:,z], Consi[:,z], elli[:,z], dis_toti[:,z], yi[:,z], Ai[:,z], Savi[:,z], gami[:,z], ri[:,z], tauCi[:,z], wi[:,z], tauWi[:,z], thetai[:,z], notreti[:,z], ppi[:,z], tauli[:,z], ivi[:,z], abi[:,z], nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi, HH_nonconvi[:,z])
    else
      lambdai[:,z], pci[:,z], Consi[:,z], elli[:,z], dis_toti[:,z], yi[:,z], Ai[:,z], Savi[:,z], HH_nonconvi[:,z] = HH(fagi, 30, 1e-10, 0.1,lambdai[:,z], pci[:,z], Consi[:,z], elli[:,z], dis_toti[:,z], yi[:,z], Ai[:,z], Savi[:,z], gami[:,z], ri[:,z], tauCi[:,z], wi[:,z], tauWi[:,z], thetai[:,z], notreti[:,z], ppi[:,z], tauli[:,z], ivi[:,z], abi[:,z], nagi, rhoi, sigmai, parlv0i, parlv1i, sigLi, HH_nonconvi[:,z])
    end
  end
  
  return lambdai, pci, Consi, elli, dis_toti, yi, Ai, Savi, HH_nonconvi
  
end