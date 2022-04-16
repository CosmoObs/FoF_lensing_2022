#!/usr/bin/python
# =====================================================
# Authors:
# Eduardo Valadão - eduardovaladao98@gmail.com
# =====================================================

""" Module containing the lens equations of the Point Lens, SIS (Singular Isothermal Sphere), 
SIE (Singular Isothermal Ellipsoid), SIEP (Singular Isothermal Elliptical Potential) and 
NFW (Navarro-Frenk-White) lens with optional external shear and external convergence. 
To make a gravitational lensing simulation we simply add the lens equation of a given lens 
model in the elliptical radial profile representing the source-galaxy, then we just put the 
result in a surface brightness radial profile appropriate to a source-galaxy. Furthermore, 
in the end of this module we have integrals over the areas of the SIS gravitational arcs 
related to the normalization of the surface brightness profile so that the total luminosity 
of the image is equal to the total signal. This is a work in progress, so my apologies to the 
user if any errors arise and if my code is too ugly. """

import math 
import mpmath as mp
from scipy.integrate import dblquad, quad
from scipy.special import gamma, factorial
from convert_arcsec_pix import convert_pix_arcsec

#========================================================================================================
# Lens equations and elliptical radial profile:
#========================================================================================================

def radius(x1, x2, x1_0, x2_0):
    ''' coordinate transformation from r to cartesian coordinates (x1, x2) '''
    return mp.sqrt((x1 - x1_0)**2 + (x2 - x2_0)**2) 

def elliptical_source(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' lens equation in cartesian coordinates of a ellipse with optional 
    constant external shear and external convergence '''
    y1 = (1 - k - u)*(x1 - x1_0)
    y2 = (1 - k + u)*(x2 - x2_0)
    return y1,y2

def sis(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' lens equation in cartesian coordinates of the SIS model with 
    optional constant external shear and external convergence'''
    theta = radius(x1, x2, x1_0, x2_0)
    y1 = ((1 - k - u) - (theta_e/theta))*(x1 - x1_0)
    y2 = ((1 - k + u) - (theta_e/theta))*(x2 - x2_0)
    return y1,y2

def point_lens(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' lens equation in cartesian coordinates of the Point Lens model 
    with optional constant external shear and external convergence '''
    theta = radius(x1, x2, x1_0, x2_0)
    y1 = ((1 - k - u) - (theta_e/theta)**2)*(x1 - x1_0)
    y2 = ((1 - k + u) - (theta_e/theta)**2)*(x2 - x2_0)
    return y1,y2

def sie(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' lens equation in cartesian coordinates of the SIE model 
    with optional constant external shear and external convergence'''
    fe = 1 - f
    ff = mp.sqrt(1 - fe**2)
    theta = radius(x1, x2, x1_0, x2_0)
    y1 = (1 - k - u)*(x1 - x1_0) - (theta_e)*(mp.sqrt(fe)/ff)*mp.asinh((ff*(x1 - x1_0))/(fe*theta))
    y2 = (1 - k + u)*(x2 - x2_0) - (theta_e)*(mp.sqrt(fe)/ff)*mp.asin((ff*(x2 - x2_0))/theta)
    return y1,y2

def siep(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' lens equation in cartesian coordinates of the SIEP model 
    with optional constant external shear and external convergence'''
    y1 = ((1 - k - u) - (theta_e)*mp.sqrt((1 + f)/((x1 - x1_0)**2 + ((1 - f)/(1 + f))*((x2 - x2_0)**2))))*(x1 - x1_0)
    y2 = ((1 - k + u) - (theta_e)*mp.sqrt((1 - f)/((x2 - x2_0)**2 + ((1 + f)/(1 - f))*((x1 - x1_0)**2))))*(x2 - x2_0) 
    return y1,y2 

def nfw(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' lens equation in cartesian coordinates of the NFW lens model 
    with optional constant external shear and external convergence '''
    theta = radius(x1, x2, x1_0, x2_0)
    if theta < theta_e:
        y1 = ((1 - k - u) - (4*ks)*((theta_e/theta)**2)*(mp.ln(theta/(2*theta_e)) + (2*theta_e/(mp.sqrt((theta_e**2 - theta**2))))*mp.re(mp.atanh(mp.sqrt((theta_e - theta)/(theta_e + theta))))))*(x1 - x1_0) 
        y2 = ((1 - k + u) - (4*ks)*((theta_e/theta)**2)*(mp.ln(theta/(2*theta_e)) + (2*theta_e/(mp.sqrt((theta_e**2 - theta**2))))*mp.re(mp.atanh(mp.sqrt((theta_e - theta)/(theta_e + theta))))))*(x2 - x2_0)
        return y1,y2
    elif theta == theta_e:
        y1 = ((1 - k - u) - (4*ks)*((theta_e/theta)**2)*(mp.ln(1/2) + 1))*(x1 - x1_0)
        y2 = ((1 - k - u) - (4*ks)*((theta_e/theta)**2)*(mp.ln(1/2) + 1))*(x2 - x2_0)
        return y1,y2
    elif theta > theta_e:
        y1 = ((1 - k - u) - (4*ks)*((theta_e/theta)**2)*(mp.ln(theta/(2*theta_e)) + (2*theta_e/(mp.sqrt((theta**2 - theta_e**2))))*mp.re(mp.atan(mp.sqrt((theta - theta_e)/(theta_e + theta))))))*(x1 - x1_0) 
        y2 = ((1 - k + u) - (4*ks)*((theta_e/theta)**2)*(mp.ln(theta/(2*theta_e)) + (2*theta_e/(mp.sqrt((theta**2 - theta_e**2))))*mp.re(mp.atan(mp.sqrt((theta - theta_e)/(theta_e + theta))))))*(x2 - x2_0)
        return y1,y2
    

def model_choice(model, x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' this function's output is the user chosen lens model '''
    if model == 1:
        y1,y2 = elliptical_source(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)
    elif model == 2:
        y1,y2 = sis(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)
    elif model == 3:
        y1,y2 = point_lens(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)
    elif model == 4:
        y1,y2 = sie(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)
    elif model == 5:
        y1,y2 = siep(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)
    elif model == 6:
        y1,y2 = nfw(x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)    
    return y1,y2    

def elliptical_radial_profile(model, x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks):
    ''' the elliptical radial profile is basically the famous ellipse equation but with more
    generality, since our elliptical source can have arbitrary position. This equation can be
    rewritten so that its form is basically R0^2 = (x/a)^2 + (y/b)^2, where R0 is the effective
    size of our source-galaxy. Making the coordinate transformation between (y1,y2) to (x1,x2)
    through the lens equation of a particular model gives us the desired solution to do the 
    astronomical simulation of a lensing event '''
    s11 = s1 - x1_0
    s22 = s2 - x2_0 
    y1,y2 = model_choice(model, x1, x2, s1, s2, phie, e, f, k, u, theta_e, x1_0, x2_0, ks)   
    return mp.sqrt((1 - e)*((y1 - s11)*mp.cos(phie) + (y2 - s22)*mp.sin(phie))**2 + (1 + e)*((y2 - s22)*mp.cos(phie) - (y1 - s11)*mp.sin(phie))**2)

#========================================================================================================
# Functions related to the pixelation of our model, this part has yet to be fixed:
#========================================================================================================

def pixel_size_sis(epsilon, r_e, n, b_n, e, phie, s1, s2, x1_0, x2_0):
    theta0 = math.atan2((s2 - x2_0), (s1 - x1_0))
    dx = abs(epsilon*((n*r_e)/(b_n*(1 - e*mp.cos(2*(theta0 - phie)))**(1/(2*n)))))
    return dx 

def pixel_size_point_lens(epsilon, r_e, n, b_n, e, phie, s1, s2, x1_0, x2_0):
    return pixel_size_sis(epsilon, r_e, n, b_n, e, phie, s1, s2, x1_0, x2_0)/2 

def I_0_integral_source(r_e, n, b_n, e):
    if type(n) == int:
        return (2*mp.pi*(r_e**2)*n*factorial(2*n - 1))/((b_n**(2*n))*mp.sqrt(1 - e**2))
    else:
        return (2*mp.pi*(r_e**2)*n*gamma(2*n))/((b_n**(2*n))*mp.sqrt(1 - e**2)) 

#========================================================================================================
# Functions related to the normalization of the surface brightness radial profile of the SIS model, 
# I apologize to the user again for the messy last function of this part, in the near future I hope 
# to improve it:
#========================================================================================================

def xmenos(phi, r0, s, phie, theta, e, theta_e):
    ''' inner part of the SIS + elliptical source images '''
    return theta_e + (1/(1 - e*mp.cos(2*(phi - phie))))*(s*mp.cos(theta - phi) - e*s*mp.cos(theta + phi - 2*phie) - mp.sqrt((r0**2)*(1 - e*mp.cos(2*(phi - phie))) - (s**2)*(1 - e**2)*(mp.sin(theta - phi)**2)))
def xmais(phi, r0, s, phie, theta, e, theta_e):
    ''' outer part of the SIS + elliptical source images '''
    return theta_e + (1/(1 - e*mp.cos(2*(phi - phie))))*(s*mp.cos(theta - phi) - e*s*mp.cos(theta + phi - 2*phie) + mp.sqrt((r0**2)*(1 - e*mp.cos(2*(phi - phie))) - (s**2)*(1 - e**2)*(mp.sin(theta - phi)**2)))

def sthy(r0, phie, e): 
    ''' limit where the elliptical source is tangent to the y2 axis '''
    return r0*mp.sqrt((1 + e*mp.cos(2*phie))/(1 - e**2))
def sth0(r0, phie, e):
    ''' limit where the origin is inside the elliptical source and an einstein ring arises '''
    return r0/mp.sqrt(1 - e*mp.cos(2*phie)) 

def a(e): 
    ''' semi-major axis of the elliptical source '''
    return 1/mp.sqrt(1 - e)
def b(e):
    ''' semi-minor axis of the elliptical source '''
    return 1/mp.sqrt(1 + e)

def P1(r0, s, phie, theta, e):
    ''' function related to the angle of the arcs' extremes '''
    return (a(e)**2 - b(e)**2)*(r0**2)*((a(e)**2 + b(e)**2)*(r0**2)*mp.cos(2*phie) - (s**2)*(2*mp.cos(2*(theta - phie)) + mp.cos(2*phie))) + ((a(e)**2 - b(e)**2)**2)*(r0**4) + s**4 + (s**2)*(s**2 - (a(e)**2 + b(e)**2)*(r0**2))*mp.cos(2*theta)
def P2(r0, s, phie, theta, e):
    ''' function related to the angle of the arcs' extremes '''
    return ((a(e)**2 + b(e)**2)*(s**2) - (a(e)**2 - b(e)**2)*(s**2)*mp.cos(2*(theta - phie)) - 2*((a(e)*b(e)*r0)**2))*((r0*(s**2)*mp.sin(2*theta) - (a(e)**2 - b(e)**2)*(r0**3)*mp.sin(2*phie))**2)
def P3(r0, s, phie, theta, e):
    ''' function related to the angle of the arcs' extremes '''
    return ((a(e)**2 - b(e)**2)**2)*(r0**4) + s**4 - 2*(a(e)**2 - b(e)**2)*(r0**2)*(s**2)*mp.cos(2*(theta - phie))

def phi_1(r0, s, phie, theta, e):
    ''' angle of the arcs' extremes, depending on the position and orientation of the source '''
    return mp.acos(mp.sqrt((P1(r0, s, phie, theta, e) + mp.sqrt(2*P2(r0, s, phie, theta, e)))/(2*P3(r0, s, phie, theta, e))))
def phi_2(r0, s, phie, theta, e):
    ''' angle of the arcs' extremes, depending on the position and orientation of the source '''
    return mp.acos(-mp.sqrt((P1(r0, s, phie, theta, e) - mp.sqrt(2*P2(r0, s, phie, theta, e)))/(2*P3(r0, s, phie, theta, e))))
def phi_3(r0, s, phie, theta, e):
    ''' angle of the arcs' extremes, depending on the position and orientation of the source '''
    return mp.acos(-mp.sqrt((P1(r0, s, phie, theta, e) + mp.sqrt(2*P2(r0, s, phie, theta, e)))/(2*P3(r0, s, phie, theta, e))))
def phi_4(r0, s, phie, theta, e):
    ''' angle of the arcs' extremes, depending on the position and orientation of the source '''
    return mp.acos(mp.sqrt((P1(r0, s, phie, theta, e) - mp.sqrt(2*P2(r0, s, phie, theta, e)))/(2*P3(r0, s, phie, theta, e))))

def I_0_integral_sis(r_e, n, b_n, s1, s2, phie, e, theta_e, x1_0, x2_0, r0):
    ''' calculation of the integral present in the normalization of the surface brightness 
    profile for the SIS + elliptical source model '''
    mp.dps = 100; mp.pretty = True
    s = radius(s1, s2, x1_0, x2_0) 
    theta = 0
    s0 = sth0(r0, phie, e)
    sy = sthy(r0, phie, e)
    f = lambda phi: mp.quad(lambda r: r*mp.exp(-b_n*((mp.sqrt((1 - e)*(((r - theta_e)*mp.cos(phi) - s)*mp.cos(phie) + ((r - theta_e)*mp.sin(phi))*mp.sin(phie))**2 + (1 + e)*(((r - theta_e)*mp.sin(phi))*mp.cos(phie) - ((r - theta_e)*mp.cos(phi) - s)*mp.sin(phie))**2))/r_e)**(1/n)), [xmenos(phi, r0, s, phie, theta, e, theta_e), xmais(phi, r0, s, phie, theta, e, theta_e)])
    if 0 <= phie <= mp.pi/2 or mp.pi < phie <= 3*(mp.pi/2):
        if s0 > s:
            a = mp.quad(f, [0, 2*mp.pi])
            return mp.mpmathify(a)    
        elif s0 <= s <= sy:
            a = mp.quad(f, [-phi_2(r0, s, phie, theta, e), phi_1(r0, s, phie, theta, e)])
            b = mp.quad(f, [-phi_2(r0, s, phie, theta, e) + mp.pi, phi_1(r0, s, phie, theta, e) + mp.pi])
            exx = mp.re(a)
            inn = mp.re(b)
            return mp.mpmathify(exx + inn)
        elif sy < s:
            a = mp.quad(f, [-phi_4(r0, s, phie, theta, e), phi_1(r0, s, phie, theta, e)])
            b = mp.quad(f, [-phi_4(r0, s, phie, theta, e) + mp.pi, phi_1(r0, s, phie, theta, e) + mp.pi])
            exx = mp.re(a)
            inn = mp.re(b)
            return mp.mpmathify(exx + inn) 
    elif mp.pi/2 < phie <= mp.pi or 3*(mp.pi/2) < phie <= 2*mp.pi:
        if s0 > s:
            a = mp.quad(f, [0, 2*mp.pi])
            return mp.mpmathify(a)
        elif s0 <= s <= sy:
            a = mp.quad(f, [-phi_1(r0, s, phie, theta, e), phi_2(r0, s, phie, theta, e)])
            b = mp.quad(f, [-phi_1(r0, s, phie, theta, e) + mp.pi, phi_2(r0, s, phie, theta, e) + mp.pi])
            exx = mp.re(a)
            inn = mp.re(b)
            return mp.mpmathify(exx + inn)
        elif sy < s:
            a = mp.quad(f, [-phi_1(r0, s, phie, theta, e), phi_4(r0, s, phie, theta, e)])
            b = mp.quad(f, [-phi_1(r0, s, phie, theta, e) + mp.pi, phi_4(r0, s, phie, theta, e) + mp.pi])
            exx = mp.re(a)
            inn = mp.re(b)
            return mp.mpmathify(exx + inn)

