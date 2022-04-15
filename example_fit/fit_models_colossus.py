import numpy as np
from colossus.cosmology import cosmology  
params = {'flat': True, 'H0': 70.0, 'Om0': 0.25, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.95}
cosmology.addCosmology('MICE', params)
cosmo = cosmology.setCosmology('MICE')
from colossus.halo import profile_nfw
from colossus.halo import profile_einasto


class Delta_Sigma_fit:
	# R en Mpc, Sigma M_Sun/Mpc2
	

    def __init__(self,R,DSigma,err,z,model='NFW',Min=1.e13,cin=3.,fit_alpha=True,RIN=0,ROUT=30.):

        Min = np.max([1.e11,Min])
        Min = np.min([1.e15,Min])
        cin = np.max([1,cin])
        cin = np.min([10,cin])

        xplot = R
        
        m = (R > RIN)*(R < ROUT)
        
        if model == 'NFW':

            pmodel = profile_nfw.NFWProfile(M = Min, c = cin, z = z, mdef = '200c')
            
        elif model == 'Einasto':
            
            pmodel = profile_einasto.EinastoProfile(M = Min, c = cin, z = z, mdef = '200c')
            

        try:
            
            BIN= m.sum()
            
            if model == 'NFW':
                
                out = pmodel.fit(R[m]*1000., DSigma[m]*(1.e3**2), 'DeltaSigma', q_err = err[m]*(1.e3**2), tolerance = 1.e-04,verbose=False)
                
                rhos,rs = out['x']
                prof = profile_nfw.NFWProfile(rhos = rhos, rs = rs)
                alpha = -999.
                Ndoff = float(BIN - 2)
                
            elif model == 'Einasto':
                
                if fit_alpha:
                    out = pmodel.fit(R[m]*1000., DSigma[m]*(1.e3**2), 'DeltaSigma', q_err = err[m]*(1.e3**2), tolerance = 1.e-04,verbose=False)
                    rhos,rs,alpha = out['x']
                    Ndoff = float(BIN - 3)
                else:
                    out = pmodel.fit(R[m]*1000., DSigma[m]*(1.e3**2), 'DeltaSigma', q_err = err[m]*(1.e3**2), tolerance = 1.e-04,verbose=False,mask=[True,True,False])
                    rhos,rs = out['x']
                    alpha = p.par['alpha']
                    Ndoff = float(BIN - 2)
                
                prof = profile_einasto.EinastoProfile(rhos = rhos, rs = rs, alpha = alpha)
                
            
            ajuste = prof.deltaSigma(R*1000.)/(1.e3**2)
            yplot  = prof.deltaSigma(xplot*1000.)/(1.e3**2)
            
            
            chi2=((((ajuste[m]-DSigma[m])**2)/(err[m]**2)).sum())/float(Ndoff)
            res=np.sqrt(((((np.log10(ajuste[m])-np.log10(DSigma[m]))**2)).sum())/Ndoff)

            M200 = prof.MDelta(z,'200c')
            c200 = prof.RDelta(z,'200c')/rs
        
        except:
            
            print('WARNING: the fit was not performed')
            yplot = xplot
            res   = -999.
            M200  = -999.
            c200  = -999.
            alpha = -999.
            chi2  = -999.
        

        self.xplot = xplot
        self.yplot = yplot
        self.res  = res
        self.M200 = M200
        self.c200 = c200
        self.alpha = alpha
        self.chi2 = chi2


def Delta_Sigma_NFW_2h(R,z,M200,c200,cosmo_params=params,terms='1h'):
    from colossus.lss import bias
    from colossus.halo import profile_nfw
    from colossus.halo import profile_outer
    from colossus.cosmology import cosmology  
    cosmology.addCosmology('MyCosmo', cosmo_params)
    cosmo = cosmology.setCosmology('MyCosmo')

    b = bias.haloBias(M200, model = 'tinker10', z = z, mdef = '200c')
    
    outer_term = profile_outer.OuterTermCorrelationFunction(z = z, bias = b)
    pNFW = profile_nfw.NFWProfile(M = M200, mdef = '200c', z = z, c = c200, outer_terms = [outer_term])    
    
    # Outer term integrated up to 50Mpc (Luo et al. 2017, Niemic et al 2017)
    if terms == '1h':
        ds_in  = pNFW.deltaSigmaInner(R*1.e3)
        ds = ds_in
    elif terms == '2h':
        ds_out = pNFW.deltaSigmaOuter(R*1.e3, interpolate=False, interpolate_surface_density=False, accuracy=0.01, max_r_integrate=100e3)
        ds = ds_out
    elif terms == '1h+2h':
        ds_in  = pNFW.deltaSigmaInner(R*1.e3)
        ds_out = pNFW.deltaSigmaOuter(R*1.e3, interpolate=False, interpolate_surface_density=False, accuracy=0.01, max_r_integrate=100e3)
        ds = ds_in + ds_out
    
    return ds/(1.e3**2)
