from fenics import *
import numpy, time, random, sys, os
from prettytable import PrettyTable
import gfort2py as gf
import matplotlib.pyplot as plt
from pdb import set_trace
        
def entire_bd(x, on_boundary):
    return on_boundary        

def tag_internal_boundary():
  
  tdim = mesh.topology().dim()
  mesh.init(tdim-1, tdim)                        # Creates connectivities between facets and cells
  facet_to_cell = mesh.topology()(tdim-1, tdim)  # MeshConnectivity object (stores connections between facets and cells)
  domain_values = domains.array()                # For each cell, assumes 1 or 0 (depending on subdomain)
  facet_values = int_boundary.array()            # For each facet, it will assume 1 if is an internal boundary and 0 otherwise

  
  for facet in range(len(facet_values)):
      cells = facet_to_cell(facet)                 # Returns the index of the cells which are connected to this facet
    
      if len(cells) == 2:                           # If the facet is conected with two cells, then it is an internal facet
          values = domain_values[cells]
          if values[0] != values[1]:
            facet_values[facet] = numpy.max(values) + (nsites-1)*numpy.min(values)
            
    
def tag_segment(l, v, w):
    tdim = mesh.topology().dim()
    # Creates connectivities between facets and cells
    mesh.init(tdim-1, tdim)
    # MeshConnectivity object (stores connections between facets and cells)                        
    facet_to_cell = mesh.topology()(tdim-1, tdim)  
    # For each facet, it will assume 1 if is an internal boundary and 0 otherwise
    facet_values = intboundaries.array()            
    p = numpy.zeros((2,2))
    fct = 0
    for facet in facets(mesh):
        # Returns the index of the cells which are connected to this facet
        cells = facet_to_cell(fct)                 
        k = 0
        if len(cells) == 1:
            for vertex in vertices(facet):
                p[0:2,k] = numpy.array([vertex.point().array()[0], vertex.point().array()[1]])
                k = k + 1
            
            a = min(v[1], w[1]); b = max(v[1], w[1])
            #Left
            if l == 1 and p[0,0] == 0.0 and p[0,1] == 0.0 and ((p[0:2,0]+p[0:2,1])/2)[1] >= a and ((p[0:2,0]+p[0:2,1])/2)[1] <= b:
                facet_values[fct] = int(l)

            #Right
            elif l == 2 and p[0,0] == 1.0 and p[0,1] == 1.0 and ((p[0:2,0]+p[0:2,1])/2)[1] >= a and ((p[0:2,0]+p[0:2,1])/2)[1] <= b:
                facet_values[fct] = int(l)
            

            a = min(v[0], w[0]); b = max(v[0], w[0])
            #Top
            if l == 3 and p[1,0] == 1.0 and p[1,1] == 1.0 and ((p[0:2,0]+p[0:2,1])/2)[0] >= a and ((p[0:2,0]+p[0:2,1])/2)[0] <= b:
                facet_values[fct] = int(l)

            #Botton
            elif l == 4 and p[1,0] == 0.0 and p[1,1] == 0.0 and ((p[0:2,0]+p[0:2,1])/2)[0] >= a and ((p[0:2,0]+p[0:2,1])/2)[0] <= b:
                facet_values[fct] = int(l)
        fct += 1

def omega(x):
    for cell in cells(mesh):
        p1 = numpy.array([cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1]])
        p2 = numpy.array([cell.get_vertex_coordinates()[2], cell.get_vertex_coordinates()[3]])
        p3 = numpy.array([cell.get_vertex_coordinates()[4], cell.get_vertex_coordinates()[5]])
        ic = incenter(p1, p2, p3)

        kmin = 0
        distmin = numpy.linalg.norm(ic - x[0:2])
        for k in range(1, nsites):
            dist = numpy.linalg.norm(ic - x[2*k:2*k+2])
            if dist < distmin:
                kmin = k
                distmin = dist

        domains.array()[cell.index()] = kmin

def voronoi(x, draw):

    An = 4
    Ax = numpy.array([0.0, 1.0, 1.0, 0.0, 0.0])
    Ay = numpy.array([0.0, 0.0, 1.0, 1.0, 0.0])
    Aflag = numpy.array([4, 2, 3, 1, 4])

    nsites = int(len(x)/2)
    sites = numpy.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]

    vor = []    
    if nsites >= 3:
        SHARED_LIB_NAME=f'./libvoro.so'
        MOD_FILE_NAME='voro.mod'

        user_cache_dir = './gfortcache'
        if not os.path.exists(user_cache_dir):
        # Create the directory if necessary
            os.makedirs(user_cache_dir)

        x=gf.fFort(SHARED_LIB_NAME,MOD_FILE_NAME, cache_folder=user_cache_dir)
        
        nvmax = 500
        nv = int(0)
        vx = numpy.zeros(nvmax,dtype=float)
        vy = numpy.zeros(nvmax,dtype=float)
        vflag = numpy.zeros(nvmax,dtype=int)
        sstart = numpy.zeros(nsites+1,dtype=int)
        istop = int(0)
        mydict = x.voronoi(nsites,sites.T,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop)[1]
        istop = (mydict)['istop']

        if istop == 0:
            sstart = (mydict)['sstart']
            nv = (mydict)['nv']
            vx = (mydict)['vx']
            vy = (mydict)['vy']
            vectorvflag = (mydict)['vflag']

            #Ajustando os índices das células
            for i in range(len(vectorvflag)):
                if vectorvflag[i] > 0:
                    vectorvflag[i] = vectorvflag[i] - 1

            for i in range(nsites):
                cell = []
                for k in range(sstart[i] - 1,sstart[i+1] - 1):
                    cell.append([[vx[k],vy[k]],vectorvflag[k]])
                vor.append(cell)

            if draw:
                colors = qsource.astype(int)
                x.drawvor(nsites,sites.T,colors,An,Ax,Ay,sstart,nv,vx,vy)
    
    if nsites == 2:
        istop = int(0)

        x0 = sites[0,0]; y0 = sites[0,1]; x1 = sites[1,0]; y1 = sites[1,1]
        u = y1 - y0; v = x1 - x0; w = -0.5*(y1**2 - y0**2 + x1**2 - x0**2)        

        if near(u, 0.0, DOLFIN_EPS) and near(v, 0.0, DOLFIN_EPS):
            print('In voronoi, nsites=2 and the two sites are identical.')
            istop = int(-1)
        else:
            l = 0; r = 0; b = 0; t = 0
            num = -w; den = u
            if num*den >= 0.0 and abs(num) <= abs(den):
                lam1 =  num/den; l = 1
            num = -(w+v); den = u
            if num*den >= 0.0 and abs(num) <= abs(den):
                lam2 =  num/den; r = 1
            num = -(w+u); den = v
            if num*den >= 0.0 and abs(num) <= abs(den):
                lam3 =  num/den; t = 1
            num = -w; den = v
            if num*den >= 0.0 and abs(num) <= abs(den):
                lam4 =  num/den; b = 1
            
            if (l + r + b + t) > 2:
                print('In voronoi, a reta passa por um vértice do quadrado')
                print('l, r, t e b valem:', l, r, t, b)

            if l == 1 and r == 1:
                if (v*sites[0,0] + u*sites[0,1] + w)*w > 0.0:
                    vor.append([[[0.0,0.0],-4], [[1.0,0.0],-2], [[1.0,lam2],1], [[0.0,lam1],-1]])
                    vor.append([[[0.0,lam1],0], [[1.0,lam2],-2], [[1.0,1.0],-3], [[0.0,1.0],-1]])
                else:
                    vor.append([[[0.0,lam1],1], [[1.0,lam2],-2], [[1.0,1.0],-3], [[0.0,1.0],-1]])
                    vor.append([[[0.0,0.0],-4], [[1.0,0.0],-2], [[1.0,lam2],0], [[0.0,lam1],-1]])
            elif t == 1 and b == 1:
                if (v*sites[0,0] + u*sites[0,1] + w)*w > 0.0:
                    vor.append([[[0.0,1.0],-1], [[0.0,0.0],-4], [[lam4,0.0],1], [[lam3,1.0],-3]])
                    vor.append([[[lam4,0.0],-4], [[1.0,0.0],-2], [[1.0,1.0],-3], [[lam3,1.0],0]])
                else:
                    vor.append([[[lam4,0.0],-4], [[1.0,0.0],-2], [[1.0,1.0],-3], [[lam3,1.0],1]])
                    vor.append([[[0.0,1.0],-1], [[0.0,0.0],-4], [[lam4,0.0],0], [[lam3,1.0],-3]])
            elif l == 1 and b == 1:
                if (v*sites[0,0] + u*sites[0,1] + w)*w > 0.0:
                    vor.append([[[0.0,0.0],-4], [[lam4,0.0],1], [[0.0,lam1],-1]])
                    vor.append([[[0.0,lam1],0], [[lam4,0.0],-4], [[1.0,0.0],-2], [[1.0,1.0],-3], [[0.0,1.0],-1]])
                else:
                    vor.append([[[0.0,lam1],1], [[lam4,0.0],-4], [[1.0,0.0],-2], [[1.0,1.0],-3], [[0.0,1.0],-1]])
                    vor.append([[[0.0,0.0],-4], [[lam4,0.0],0], [[0.0,lam1],-1]])
            elif r == 1 and b == 1:
                if (v*sites[0,0] + u*sites[0,1] + w)*w > 0.0:
                    vor.append([[[0.0,0.0],-4], [[lam4,0.0],1], [[1.0,lam2],-2], [[1.0,1.0],-3], [[0.0,1.0],-1]])
                    vor.append([[[lam4,0.0],-4], [[1.0,0.0],-2], [[1.0,lam2],0]])
                else:
                    vor.append([[[lam4,0.0],-4], [[1.0,0.0],-2], [[1.0,lam2],1]])
                    vor.append([[[0.0,0.0],-4], [[lam4,0.0],0], [[1.0,lam2],-2], [[1.0,1.0],-3], [[0.0,1.0],-1]])
            elif r == 1 and t == 1:
                if (v*sites[0,0] + u*sites[0,1] + w)*w > 0.0:
                    vor.append([[[0.0,0.0],-4], [[1.0,0.0],-2], [[1.0,lam2],1], [[lam3,1.0],-3], [[0.0,1.0],-1]])
                    vor.append([[[1.0,lam2],-2], [[1.0,1.0],-3], [[lam3,1.0],0]])
                else:
                    vor.append([[[1.0,lam2],-2], [[1.0,1.0],-3], [[lam3,1.0],1]])
                    vor.append([[[0.0,0.0],-4], [[1.0,0.0],-2], [[1.0,lam2],0], [[lam3,1.0],-3], [[0.0,1.0],-1]])
            elif l == 1 and t == 1:
                if (v*sites[0,0] + u*sites[0,1] + w)*w > 0.0:
                    vor.append([[[0.0,0.0],-4], [[1.0,0.0],-2], [[1.0,1.0],-3], [[lam3,1.0],1], [[0.0,lam1],-1]])
                    vor.append([[[0.0,lam1],0], [[lam3,1.0],-3], [[0.0,1.0],-1]])
                else:
                    vor.append([[[0.0,lam1],1], [[lam3,1.0],-3], [[0.0,1.0],-1]])
                    vor.append([[[0.0,0.0],-4], [[1.0,0.0],-2], [[1.0,1.0],-3], [[lam3,1.0],0], [[0.0,lam1],-1]])
        
    return vor, istop

def rotate(x):
    return numpy.array([-x[1], x[0]])

def psiv(a, v, w):
    return Expression("-((x[0] - w0)*aw0 + (x[1] - w1)*aw1) / denpsiv", degree=2, w0 = w[0], w1 = w[1], aw0 = (rotate(a - w))[0], aw1 = (rotate(a - w))[1], denpsiv = numpy.inner(a - v, rotate(a - w)))

def gradpsiv(a, v, w):
    den = numpy.inner(a - v, rotate(a - w))
    if near(den, 0.0, DOLFIN_EPS):
        print('In gradpsiv, den is close to zero.')
    return - rotate(a - w) / den

def Mvl(v, gradphil, aj, ai):
    detval = numpy.linalg.det(numpy.array([aj - ai, gradphil]))
    if near(detval, 0.0, DOLFIN_EPS):
        print('In Mvl, detval is close to zero')
    return - numpy.outer(rotate(gradphil), v - ai) / detval

def Mv(v,ai, aj, ak):
    detval = numpy.linalg.det(numpy.array([aj - ai, ak - ai]))
    if near(detval, 0.0, DOLFIN_EPS):
        print('In Mv, detval is close to zero')
    return numpy.outer(rotate(ai - aj), (v - ak)) / detval

def tcoeff(p0, p1, p2):
    tcoeffval = numpy.zeros((3,3))

    T = numpy.array([p0, p1, p2])
    
    for i in range(3):
        ip1 = (i + 1) % 3
        ip2 = (i + 2) % 3
    
        a = T[ip1,1] - T[i,1]
        b = T[i,0] - T[ip1,0]
        c = (T[i,1] - T[ip1,1]) * T[i,0] + (T[ip1,0] - T[i,0]) * T[i,1]

        tcoeffval[i,:] = [a,b,c]

        if a * T[ip2,0] + b * T[ip2,1] + c <= 0.0:
            tcoeffval[i,:] = - tcoeffval[i,:]

    return tcoeffval

def triangintdom(site, v, vprev, vnext):
    # # Coefficients of the straight lines that determine a triangle.
    tcoeffvalprev = tcoeff(site,v,vprev)
    a0 = tcoeffvalprev[0,0]; b0 = tcoeffvalprev[0,1]; c0 = tcoeffvalprev[0,2]
    a1 = tcoeffvalprev[1,0]; b1 = tcoeffvalprev[1,1]; c1 = tcoeffvalprev[1,2]
    a2 = tcoeffvalprev[2,0]; b2 = tcoeffvalprev[2,1]; c2 = tcoeffvalprev[2,2]

    tcoeffvalnext = tcoeff(site,v,vnext)
    A0 = tcoeffvalnext[0,0]; B0 = tcoeffvalnext[0,1]; C0 = tcoeffvalnext[0,2]
    A1 = tcoeffvalnext[1,0]; B1 = tcoeffvalnext[1,1]; C1 = tcoeffvalnext[1,2]
    A2 = tcoeffvalnext[2,0]; B2 = tcoeffvalnext[2,1]; C2 = tcoeffvalnext[2,2]

    for cell in cells(mesh):
        p1 = numpy.array([cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1]])
        p2 = numpy.array([cell.get_vertex_coordinates()[2], cell.get_vertex_coordinates()[3]])
        p3 = numpy.array([cell.get_vertex_coordinates()[4], cell.get_vertex_coordinates()[5]])
        ic = incenter(p1, p2, p3)

        if a0 * ic[0] + b0 * ic[1] + c0 >= 0.0 and a1 * ic[0] + b1 * ic[1] + c1 >= 0.0 and a2 * ic[0] + b2 * ic[1] + c2 >= 0.0:
            triangdomains.array()[cell.index()] = 1
        elif A0 * ic[0] + B0 * ic[1] + C0 >= 0.0 and A1 * ic[0] + B1 * ic[1] + C1 >= 0.0 and A2 * ic[0] + B2 * ic[1] + C2 >= 0.0:
            triangdomains.array()[cell.index()] = 2
        else: 
            triangdomains.array()[cell.index()] = 0

def incenter(A, B, C):
    a = numpy.linalg.norm(B-C)
    b = numpy.linalg.norm(C-A)
    c = numpy.linalg.norm(A-B)
    return numpy.array([(a*A[0]+b*B[0]+c*C[0])/(a+b+c), (a*A[1]+b*B[1]+c*C[1])/(a+b+c)])

def xinit(ninit):
    for j in range(ninit):
        # Perturbando a solução e projetando em D = [0,1]x[0,1]
        for i in range(2*nsites):
            xini[i] = max(0.0, min(1.0, solx[i] + 0.1*(2*random.random() - 1)))
    return xini

def noise():
    noise_num = 0.0 
    noise_denom = 0.0

    for alpha in range(nsources):
        zeta_alpha = 1.0 / weightsG[alpha]
        hsolve_clean = Function(V1)
        hsolve = Function(V1)
        htrial = TrialFunction(V1)
        vtest = TestFunction(V1)
        # Define PDE
        a = inner(grad(htrial), grad(vtest)) * dx
        for i in range(nsites):
            a = a + Constant(qsource[i])*htrial * vtest * dx(i)     
        L = fsource[alpha] * vtest * dx                 
        # Assemble matrices
        aa = assemble(a)
        LL = assemble(L)
        # Dirichlet boundary conditions  
        DirichletBC(V1,0.0,entire_bd).apply(aa)
        DirichletBC(V1,0.0,entire_bd).apply(LL)       
        solve(aa, hsolve_clean.vector(), LL)  
        # add noise to the data
        max_hsolve = numpy.abs(hsolve_clean.vector()[:]).max()
        h_perturb = numpy.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
        hsolve.vector()[:] = hsolve_clean.vector()[:] + h_perturb
        noise_num = noise_num + zeta_alpha*assemble( (hsolve - hsolve_clean)**2 * dx ) 
        noise_denom = noise_denom +  zeta_alpha*assemble( hsolve_clean**2 * dx )
    
    noise_level = numpy.sqrt(noise_num / noise_denom)
    return noise_level

def evalfg(n, x, vor, greq=False, gvol=False, gbound=False):
    nsites = int(n/2)
    sites = numpy.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]
    
    omega(x)
    dx = Measure("dx", domain=mesh, subdomain_data=domains)

    int_boundary.set_all(0)
    tag_internal_boundary()
    dS_int = Measure("dS", domain=mesh, subdomain_data=int_boundary)

    G = 0.0
    gradGvol = numpy.zeros(2*nsites)
    gradGbound = numpy.zeros(2*nsites)

    for alpha in range(nsources):
 
        # Define rhs and lhs of PDE (8)
        a = inner(grad(utf), grad(wtest)) * dx
        for i in range(nsites):     
            a = a + Constant(qsource[i])*utf * wtest *  dx(i)
        L = fsource[alpha] * wtest * dx          
        # Assemble matrices
        aa = assemble(a)
        LL = assemble(L)
        # Dirichlet boundary conditions
        DirichletBC(V1,0.0,entire_bd).apply(aa)
        DirichletBC(V1,0.0,entire_bd).apply(LL)                 
        # Solve
        solve(aa, usol.vector(), LL)

        # Compute cost function        
        Galpha = assemble( (usol - h[alpha])**2 * dx )
        if weightsG[alpha] == -1.0:
            weightsG[alpha] = 1.0/Galpha
    
        G = G + 0.5 * weightsG[alpha] * Galpha
        # print('Galpha =', Galpha)
        
        if greq:
            # Define rhs and lhs of PDE (?)
            a = inner(grad(ptf), grad(wtest)) * dx
            for i in range(nsites):     
                a = a + Constant(qsource[i])*ptf * wtest *  dx(i)
            L = (h[alpha] - usol) * wtest * dx  
            # Assemble matrices
            aa = assemble(a)
            LL = assemble(L)
            # Dirichlet boundary conditions        
            DirichletBC(V1,0.0,entire_bd).apply(aa)
            DirichletBC(V1,0.0,entire_bd).apply(LL)   
            solve(aa, psol.vector(), LL)

            for i in range(nsites):
                if gvol:
                    ### Volumetric Approach ###
                    temp = 0.5*(usol - h[alpha])**2 - fsource[alpha]*psol + inner(grad(usol), grad(psol)) + Constant(qsource[i])*usol*psol

                    # Warning: In case of non-constant fsource and Volumetric Approach, the derivative of fsource must be calculated explicitly.
                    # S0 = (h[alpha] - usol)*grad(h[alpha]) - grad(fsource[alpha])
                    S0 = (h[alpha] - usol)*grad(h[alpha])

                    S1 = Identity(2)*temp
                    S1 = S1 - outer(grad(psol),grad(usol)) - outer(grad(usol),grad(psol))
                    
                if gbound:
                    ### Boundary Expression ###
                    S1bound = Constant(qsource[i])*usol*psol

                vprev = ((vor[i])[len(vor[i]) - 1])[0]
                vprevflag = ((vor[i])[len(vor[i]) - 1])[1]

                for r in range(len(vor[i])):
                    v = ((vor[i])[r])[0]
                    vflag = ((vor[i])[r])[1]
                    vnext = ((vor[i])[(r+1)%len(vor[i])])[0]
                    vnextflag = ((vor[i])[(r+1)%len(vor[i])])[1]

                    ### Boundary Expression ###
                    if gbound and (vflag >= 0):
                        # As células vizinhas da aresta interna E são a_i e a_vflag
                        den = numpy.linalg.norm(sites[vflag,:] - sites[i, :])
                        xi = int(numpy.max([vflag, i]) + (nsites-1)*numpy.min([vflag, i]))

                        if i > vflag:
                            hkE0 = assemble(S1bound*Expression("x[0] - ak", ak = sites[vflag,0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hkE1 = assemble(S1bound*Expression("x[1] - ak", ak = sites[vflag,1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE0 = assemble(S1bound*Expression("x[0] - ai", ai = sites[i, 0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE1 = assemble(S1bound*Expression("x[1] - ai", ai = sites[i, 1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                        else:
                            hkE0 = assemble(S1bound*Expression("x[0] - ak", ak = sites[vflag,0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hkE1 = assemble(S1bound*Expression("x[1] - ak", ak = sites[vflag,1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE0 = assemble(S1bound*Expression("x[0] - ai", ai = sites[i, 0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE1 = assemble(S1bound*Expression("x[1] - ai", ai = sites[i, 1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            
                        gradGbound[2*vflag:2*vflag+2] = gradGbound[2*vflag:2*vflag+2] - weightsG[alpha] * numpy.array([hkE0, hkE1])

                        gradGbound[2*i:2*i+2] = gradGbound[2*i:2*i+2] + weightsG[alpha] * numpy.array([hiE0, hiE1]) 
                                            
                    ### Volumetric Approach ###
                    if gvol and (vprevflag >= 0 or vflag >= 0):    
                        # Tag the integration domain
                        triangintdom(sites[i,:],v,vprev, vnext)
                        dt = Measure('dx')(subdomain_data=triangdomains)                  
                        
                        psivvalvprev = psiv(sites[i,:], v, vprev)
                        gradpsivvalTemp = gradpsiv(sites[i,:],v,vprev)
                        gradpsivval = Constant((gradpsivvalTemp[0], gradpsivvalTemp[1]))
                  
                        Tivw0 = assemble(((S1*gradpsivval + S0*psivvalvprev)[0])*dt(1))
                        Tivw1 = assemble(((S1*gradpsivval + S0*psivvalvprev)[1])*dt(1))
                        
                        psivvalvnext = psiv(sites[i,:], v, vnext)
                        gradpsivvalTemp = gradpsiv(sites[i,:],v,vnext)
                        gradpsivval = Constant((gradpsivvalTemp[0], gradpsivvalTemp[1]))                  
                    
                        Tivw0 = Tivw0 + assemble(((S1*gradpsivval + S0*psivvalvnext)[0])*dt(2))
                        Tivw1 = Tivw1 + assemble(((S1*gradpsivval + S0*psivvalvnext)[1])*dt(2))

                        if (vprevflag < 0 or vflag < 0):
                            j = max(vprevflag, vflag)
                            l = int(- min(vprevflag, vflag))

                            # Definindo qual vizinho está no intervalo de integração
                            if vflag >= 0:
                                w = vprev
                            else:
                                w = vnext

                            gradpsivvalTemp = gradpsiv(sites[i,:],v,w)
                            gradpsivval = Constant((gradpsivvalTemp[0], gradpsivvalTemp[1])) 

                            mvl = Mvl(v,gradphi[l],sites[j,:],sites[i,:])
                            gradGvol[2*i:2*i+2] = gradGvol[2*i:2*i+2] + weightsG[alpha] * numpy.matmul(mvl.T,[Tivw0, Tivw1])
                            
                            mvl = Mvl(v,gradphi[l],sites[i,:],sites[j,:])
                            gradGvol[2*j:2*j+2] = gradGvol[2*j:2*j+2] + weightsG[alpha] * numpy.matmul(mvl.T,[Tivw0, Tivw1])
                        else:
                            j = vprevflag; k = vflag #it could also be j = vflag; k = vprevflag 

                            mvjki = Mv(v,sites[j,:],sites[k,:], sites[i,:])
                            gradGvol[2*i:2*i+2] = gradGvol[2*i:2*i+2] + weightsG[alpha] * numpy.matmul(mvjki.T,[Tivw0, Tivw1])
                            
                            mvkij = Mv(v,sites[k,:],sites[i,:], sites[j,:])
                            gradGvol[2*j:2*j+2] = gradGvol[2*j:2*j+2] + weightsG[alpha] * numpy.matmul(mvkij.T,[Tivw0, Tivw1])

                            mvijk = Mv(v,sites[i,:],sites[j,:], sites[k,:])
                            gradGvol[2*k:2*k+2] = gradGvol[2*k:2*k+2] + weightsG[alpha] * numpy.matmul(mvijk.T,[Tivw0, Tivw1])
                            
                    vprev = v
                    vprevflag = vflag
                             
    if greq:
        if gvol == True and gbound == True:
            return G, gradGvol, gradGbound
        if gvol == True and gbound == False:
            return G, gradGvol
        if gvol == False and gbound == True:
            return G, gradGbound
    else:
        return G

def projectintoA(p):
    p[0] = max(0.01, min(p[0], 0.99))
    p[1] = max(0.01, min(p[1], 0.99))

def project(n, x):
    for i in range(int(n/2)):
        projectintoA(x[2*i:2*i+2])

def projectedgradient(n, x, epsg, maxit):
    project(n, x)

    vor, vorflag = voronoi(x, draw = True)
    if vorflag != 0:
        print('In projectedgradient, voronoi diagram is not well defined at (projection of) the initial guess')
        sys.exit()

    global weightsG
    weightsG = -numpy.ones(nsources)

    f, g = evalfg(n, x, vor, greq = True, gvol=False, gbound=True)   

    numevalf = 1

    gp = x - g
    project(n, gp)
    gp = gp - x
    normgp = numpy.linalg.norm(gp)

    # Saving the initial values of G(x^0) and |grad G(x^0)|
    finit, normgpinit = f,normgp

    iter = 0

    myTable = PrettyTable(["iter", "fcnt", "G", "||gP||", "x"])
    myTable.add_row([iter, numevalf, f, normgp, x])
    data = myTable.get_string()
    with open('./saida.txt', 'w') as txt:
        txt.write(data)
    print(myTable)
    smallstep = False
    while normgp > epsg and iter < maxit and not smallstep:
        iter = iter + 1

        alpha = 1.0
        xtrial = x + alpha * gp
        vor, vorflag = voronoi(xtrial, draw = False)
        if vorflag != 0:
            ftrial = float('inf')
        else:
            ftrial = evalfg(n, xtrial, vor)
            numevalf = numevalf + 1

        print('-->', alpha, ftrial, numevalf)
        
        gtgp = numpy.inner(g, gp)

        while not (ftrial <= f + 1E-4 * alpha * gtgp) and not smallstep:
            atrial = (- gtgp * alpha**2) / (2.0*(ftrial - f - alpha * gtgp))
            if not (0.1*alpha <= atrial and atrial <= 0.9*alpha):
                atrial = alpha / 2
            alpha = atrial
            xtrial = x + alpha * gp
            vor, vorflag = voronoi(xtrial, draw = False)
            if vorflag != 0:
                ftrial = float('inf')
            else:
                ftrial = evalfg(n, xtrial, vor)
                numevalf = numevalf + 1

            print('-->', alpha, ftrial, numevalf)
        
            if alpha < 1E-6:
                smallstep = True

        x = xtrial
        vor, vorflag = voronoi(x, draw = True)
        if vorflag != 0:
            print('In projectedgradient, vorflag must be zero here and it is not.')

        f, g = evalfg(n, x, vor, greq = True, gvol = False, gbound = True)
        numevalf = numevalf + 1

        gp = x - g
        project(n, gp)
        gp = gp - x
        normgp = numpy.linalg.norm(gp)

        myTable.add_row([iter, numevalf, f, normgp, x])
        data = myTable.get_string()
        with open('./saida.txt', 'w') as txt:
            txt.write(data)
        print(myTable)
             

    if iter > maxit-1:
        print('Maximum number of iterations reached')
        flagsol = 0
    elif normgp <= epsg:
        print('Small search direction')
        flagsol = 1
    elif smallstep:
        print('Too small step in line search')
        flagsol = 2
    else:
        print('In projectedgradient, main loop ended by an unknown criterion.')
    
    return flagsol, x, finit, normgpinit, f, normgp, iter, numevalf
    
    
def randomInit(ntrials, nsites):
    Gtrial = numpy.zeros(ntrials)
    xtrial_mat = numpy.random.rand(2*nsites, ntrials)  
    
    for k in range(ntrials):
        xtrial = xtrial_mat[:,k]
        omega(xtrial)
        dx = Measure("dx", domain=mesh, subdomain_data=domains)
        
        Gtrial[k] = 0.0
        for alpha in range(nsources):
            a = inner(grad(utf), grad(wtest)) * dx
            for i in range(nsites):     
                a = a + Constant(qsource[i])*utf * wtest *  dx(i)
            L = fsource[alpha] * wtest * dx          
            aa = assemble(a)
            LL = assemble(L)
            DirichletBC(V1,0.0,entire_bd).apply(aa)
            DirichletBC(V1,0.0,entire_bd).apply(LL)                 
            solve(aa, usol.vector(), LL)
        
        
            # Compute cost function        
            Galpha = assemble( (usol - h[alpha])**2 * dx )
            Gtrial[k] = Gtrial[k] + 0.5 * Galpha
    
    print('---------------------')        
    print('Random initialization')    
    print('Gtrial = ', Gtrial)
    print('---------------------')        

    kmin = numpy.argmin(Gtrial)        
    xiniRand = xtrial_mat[:,kmin]
    
    return xiniRand


def vorDiag(qsource, m, mesh, V, sites):
    # Create a MeshFunction to store values on cells
    cell_values = MeshFunction('double', mesh, dim=2)

    domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    domain_values = domains.array() 

    # Assign values to each cell
    for cell in cells(mesh):
        # The index() function returns the index of the cell
        #cell_values[cell] = cell.index()
        vertex_coordinates = cell.get_vertex_coordinates()
        # Print the coordinates
        #print(f"Coordinates of cell {cell.index()}: {vertex_coordinates}")    
    
        midpoint_x = (vertex_coordinates[0] + vertex_coordinates[2] + vertex_coordinates[4])/3.0
        midpoint_y = (vertex_coordinates[1] + vertex_coordinates[3] + vertex_coordinates[5])/3.0
       
        #print('vertex_coordinates = ',vertex_coordinates)
        dist = numpy.zeros(m) 
        for k in range(m):
            dist[k] = (midpoint_x - sites[2*k])**2 + (midpoint_y - sites[2*k+1])**2  
    
        cell_values[cell] = qsource[numpy.argmin(dist)]
        domain_values[cell.index()] = numpy.argmin(dist) # mark the subdomains with the index of the Voronoi cell
    
    f = Function(V)

    dm = V.dofmap()
    for cell in cells(mesh):
        f.vector()[dm.cell_dofs(cell.index())] = cell_values[cell] 
        
    return f


def plotVor(qsource, nsites, mesh, V, solx, xini, xcurrent, nsources, nmesh):

    input_files = './potential_loop/'+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'/'+str(ninit)+'/'
    if not os.path.exists(input_files):
        # Create the directory if necessary
        os.makedirs(input_files)

    ground = vorDiag(qsource, nsites, mesh, V, solx)  
    fig = plt.figure()
    fig.colorbar(plot(ground))
    #plt.title('Ground truth')
    plt.savefig(input_files + 'ground.png')
  
    fini = vorDiag(qsource, nsites, mesh, V, xini)  
    fig = plt.figure()
    fig.colorbar(plot(fini))
    #plt.title('Initialization')
    #plot(mesh)
    plt.savefig(input_files + 'init.png')

    fopt = vorDiag(qsource, nsites, mesh, V, xcurrent)  
    fig = plt.figure()
    fig.colorbar(plot(fopt))
    #plt.title('Optimization')
    plt.savefig(input_files + 'optimized.png')

    groundnorm = assemble( ground * dx )
    error_init = assemble( abs(fini - ground) * dx )/groundnorm
    error_opt = assemble( abs(fopt - ground) * dx )/groundnorm


    print('error_init = ', error_init)
    print('error_opt = ', error_opt)

    return error_opt, error_init 

def get_data(nsites):
    #----------------
    # Ground truth data
    #----------------
    if nsites == 2:
        solx = numpy.array([0.50, 0.25,
                            0.50, 0.75])
    if nsites == 3:
        solx = numpy.array([0.23, 0.65,
                            0.87, 0.78,
                            0.72, 0.32])
    if nsites == 4:
        solx = numpy.array([0.17, 0.89,
                            0.25, 0.20,
                            0.64, 0.41,
                            0.73, 0.93])
    if nsites == 5:
        solx = numpy.array([0.23, 0.91,
                            0.24, 0.33,
                            0.93, 0.38,
                            0.90, 0.85,
                            0.45, 0.54])
    if nsites == 6:
        solx = numpy.array([0.56, 0.67,
                            0.42, 0.37,
                            0.91, 0.10,
                            0.15, 0.50,
                            0.90, 0.86,
                            0.18, 0.90])
    if nsites == 7:
        solx = numpy.array([0.15, 0.94,
                            0.21, 0.47,
                            0.83, 0.39,
                            0.14, 0.10,
                            0.89, 0.78,
                            0.49, 0.76,
                            0.52, 0.30])
    if nsites == 8:
        solx = numpy.array([0.41, 0.29,
                            0.16, 0.54,
                            0.80, 0.48,
                            0.14, 0.10,
                            0.90, 0.85,
                            0.51, 0.77,
                            0.91, 0.21,
                            0.13, 0.96])
    if nsites == 9:
        solx = numpy.array([0.50, 0.21,
                            0.60, 0.25,
                            0.37, 0.11,
                            0.37, 0.87,
                            0.88, 0.69,
                            0.10, 0.85,
                            0.79, 0.71,
                            0.24, 0.17,
                            0.43, 0.80])
    if nsites == 10:
        solx = numpy.array([0.38, 0.78,
                            0.55, 0.33,
                            0.31, 0.84,
                            0.50, 0.49,
                            0.16, 0.36,
                            0.44, 0.15,
                            0.28, 0.69,
                            0.57, 0.13,
                            0.45, 0.44,
                            0.51, 0.89])
    if nsites == 11:    
        solx = numpy.array([0.15, 0.30,
                            0.51, 0.31,
                            0.71, 0.80,
                            0.32, 0.81,
                            0.70, 0.17,
                            0.89, 0.73,
                            0.52, 0.57,
                            0.74, 0.48,
                            0.76, 0.12,
                            0.61, 0.36,
                            0.20, 0.75])
    if nsites == 12:
        solx = numpy.array([0.89, 0.36,
                            0.70, 0.70,
                            0.31, 0.82,
                            0.16, 0.63,
                            0.50, 0.37,
                            0.71, 0.89,
                            0.42, 0.80,
                            0.87, 0.19,
                            0.25, 0.16,
                            0.60, 0.21,
                            0.27, 0.89,
                            0.35, 0.14])
    if nsites == 13:
        solx = numpy.array([0.71, 0.80, 
                            0.32, 0.81,
                            0.70, 0.17,
                            0.74, 0.69,
                            0.52, 0.57,
                            0.74, 0.48,
                            0.76, 0.12,
                            0.61, 0.36,
                            0.18, 0.73,
                            0.75, 0.37,
                            0.19, 0.34,
                            0.18, 0.19,
                            0.38, 0.86])
    if nsites == 14:
        solx = numpy.array([0.20, 0.39,
                            0.12, 0.35,
                            0.10, 0.56,
                            0.42, 0.71,
                            0.65, 0.13,
                            0.29, 0.77,
                            0.21, 0.15,
                            0.31, 0.55,
                            0.82, 0.45,
                            0.58, 0.36,
                            0.32, 0.88,
                            0.09, 0.73,
                            0.75, 0.83,
                            0.54, 0.58])
    if nsites == 15:
        solx = numpy.array([0.41, 0.29,
                            0.16, 0.54,
                            0.80, 0.48,
                            0.14, 0.10,
                            0.90, 0.85,
                            0.51, 0.80,
                            0.91, 0.21,
                            0.15, 0.90,
                            0.53, 0.06,
                            0.50, 0.50,
                            0.20, 0.75,
                            0.68, 0.97,
                            0.13, 0.34,
                            0.67, 0.27,
                            0.70, 0.07])

    if nsites == 16:
        solx = numpy.array([0.90, 0.90,
                            0.15, 0.90,
                            0.20, 0.29,
                            0.10, 0.79,
                            0.58, 0.11,
                            0.22, 0.77,
                            0.32, 0.64,
                            0.55, 0.61,
                            0.43, 0.27,
                            0.11, 0.56,
                            0.75, 0.40,
                            0.90, 0.11,
                            0.37, 0.26,
                            0.46, 0.75,
                            0.33, 0.45,
                            0.65, 0.65])

    if nsites == 17:
        solx = numpy.array([0.61, 0.85,
                            0.50, 0.10,
                            0.88, 0.50,
                            0.22, 0.86,
                            0.82, 0.62,
                            0.57, 0.50,
                            0.76, 0.53,
                            0.10, 0.35,
                            0.15, 0.16,
                            0.77, 0.40,
                            0.85, 0.20,
                            0.38, 0.45,
                            0.44, 0.84,
                            0.85, 0.80,
                            0.39, 0.70,
                            0.28, 0.59,
                            0.63, 0.73])

    if nsites == 18:
        solx = numpy.array([0.59, 0.28,
                            0.70, 0.82,
                            0.19, 0.60,
                            0.60, 0.53,
                            0.67, 0.32,
                            0.39, 0.66,
                            0.16, 0.37,
                            0.54, 0.65,
                            0.17, 0.86,
                            0.38, 0.20,
                            0.44, 0.86,
                            0.53, 0.56,
                            0.75, 0.36,
                            0.77, 0.59,
                            0.23, 0.53,
                            0.18, 0.56,
                            0.39, 0.87,
                            0.88, 0.15])

    if nsites == 19:
        solx = numpy.array([0.03, 0.42,
                            0.37, 0.03,
                            0.32, 0.90,
                            0.22, 0.58,
                            0.88, 0.14,
                            0.65, 0.83,
                            0.90, 0.54,
                            0.50, 0.70,
                            0.38, 0.23,
                            0.10, 0.75,
                            0.40, 0.50,
                            0.79, 0.96,
                            0.77, 0.56,
                            0.85, 0.23,
                            0.25, 0.45,
                            0.86, 0.46,
                            0.55, 0.23,
                            0.66, 0.46,
                            0.90, 0.74])

    return solx


    
################## 
################## 
# MAIN ALGORITHM
##################
##################

# Start Time
starttime = time.time()

nsites = int(sys.argv[1])
ninit = int(sys.argv[2])
nsources = int(sys.argv[3])
nmesh = int(sys.argv[4])
noise_coeff = float(sys.argv[5])
typeProblem = int(sys.argv[6])
typeinit = int(sys.argv[7])
Num = int(sys.argv[8])

if typeProblem == 2:
    qbinary = True
    ternary = False
elif typeProblem == 3:
    qbinary = False
    ternary = True
else:
    qbinary = False
    ternary = False


# Create unit square mesh
mesh = UnitSquareMesh(nmesh,nmesh, 'crossed')
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V1 = FunctionSpace(mesh, P1)

# Define domains
domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)

# Definindo o domínio triangular de integração
triangdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)

# Define internal interface domain
int_boundary = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

# Definindo o domínio linear de integração
intboundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

#Funções phi_l que definem o domínio D: left-> phi_1 = -x[0], right-> phi_2 = x[0] - 1, top-> phi_3 = -x[1], bottom-> phi_4 = x[1] - 1 
#Gradientes:

gradphi_1 = [-1.0, 0.0]
gradphi_2 = [1.0, 0.0]
gradphi_3 = [0.0, -1.0]
gradphi_4 = [0.0, 1.0]
gradphi = [[0.0, 0.0], gradphi_1, gradphi_2, gradphi_3, gradphi_4]

#----------------
# Data
#----------------
qsource = 10.0*numpy.arange(1.0, nsites+1)
if qbinary: 
    qsource = numpy.zeros(nsites)
    qbasemin = numpy.arange(1.0, 3.0, 0.1)
    qbasemax = numpy.arange(9.0, 11.0, 0.1)
    qsource[0:int(nsites/2)+1] = 10.0*qbasemin[0:int(nsites/2)+1]
    qsource[-int(nsites/2):] = 10.0*qbasemax[-int(nsites/2):]
if ternary:
    for l in range(nsites):
        for i in range(3):
            if l % 3 == i:
                qsource[l] = 5.0 + 5*i

fsource_list = [Expression('1.0', degree = 2), Expression("cos(pi*x[0])*cos(pi*x[1])", degree=2), Expression("sin(pi*x[0])*sin(pi*x[1])", degree=2), Expression("cos(2*pi*x[0])*cos(2*pi*x[1])", degree=2)]

if nsources == 1:
    fsource = fsource_list[0:1]
if nsources == 2:
    fsource = fsource_list[0:2]
if nsources == 3:
    fsource = fsource_list[0:3]
if nsources == 4:
    fsource = fsource_list[0:4] 

solx = get_data(nsites)

# Determinando os vértices das células de Voronoi
vor, istop = voronoi(solx, draw = False)
if istop != 0:
    print('The Voronoi method encountered an error while constructing the manufactured solution')
    sys.exit()
    
#  Tag subdomains
omega(solx)

dx = Measure("dx", domain=mesh, subdomain_data=domains)

h = []
for alpha in range(nsources):

    hsolve_clean = Function(V1)
    hsolve = Function(V1)
    htrial = TrialFunction(V1)
    vtest = TestFunction(V1)
    # Define PDE
    a = inner(grad(htrial), grad(vtest)) * dx
    for i in range(nsites):
        a = a + Constant(qsource[i])*htrial * vtest * dx(i)     
    L = fsource[alpha] * vtest * dx                 
    # Assemble matrices
    aa = assemble(a)
    LL = assemble(L)
    # Dirichlet boundary conditions  
    DirichletBC(V1,0.0,entire_bd).apply(aa)
    DirichletBC(V1,0.0,entire_bd).apply(LL)       
    #solve(a == L, hsolve)
    solve(aa, hsolve_clean.vector(), LL)  
    # add noise to the data
    max_hsolve = numpy.abs(hsolve_clean.vector()[:]).max()
    h_perturb = numpy.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
    hsolve.vector()[:] = hsolve_clean.vector()[:] + h_perturb

    # save data
    h.append(hsolve)


# Define trial functions, test functions and solutions
utf = TrialFunction(V1)
ptf = TrialFunction(V1)
wtest = TestFunction(V1)
usol = Function(V1)
psol = Function(V1)

#Initial parameter for projected gradient
maxit = 500
epsg = 1E-6
random.seed(123456)
xini = numpy.zeros(2*nsites)

if typeinit == 1:
    xini = xinit(ninit)
elif typeinit == 2:
    global weightsG
    weightsG = numpy.ones(nsources)
    xini = xinit(ninit)
    vor, istop = voronoi(xini, draw = False)
    if istop != 0:
        print('The Voronoi method encountered an error while constructing the manufactured solution')
        sys.exit()
    feval = evalfg(2*nsites, xini, vor)
    for j in range(Num - 1):
        xcurrent = xinit(ninit)
        vor, istop = voronoi(xcurrent, draw = False)
        if istop != 0:
            print('The Voronoi method encountered an error while constructing the manufactured solution')
            sys.exit()
        fevaltrial = evalfg(2*nsites, xcurrent, vor)
        if fevaltrial < feval:
            xini = numpy.copy(xcurrent)
            feval = fevaltrial
elif typeinit == 3: 
    ntrials = Num
    xini = randomInit(ntrials, nsites)
    

flagsol, xcurrent, finit, normgpinit, ffinal, normgpfinal, iter, numevalf = projectedgradient(2*nsites, xini, epsg, maxit)

## Compute noise level##
#  Tag subdomains
omega(solx)
dx = Measure("dx", domain=mesh, subdomain_data=domains) 

noise_level = noise()

# Final time
finaltime = time.time()

CPU_time = finaltime - starttime

input_files = './potential_loop/'+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'/'+str(ninit)+'/'
if not os.path.exists(input_files):
    # Create the directory if necessary
    os.makedirs(input_files)

# Save images
erroropt, errorinit = plotVor(qsource, nsites, mesh, V1, solx, xini, xcurrent, nsources, nmesh)

# Save data for drawing
numpy.savez(input_files+'drawdata.npz', qsource=qsource, solx=solx, xini=xini, xcurrent=xcurrent)
# Save table of iterations
os.rename('./saida.txt', input_files + 'saida.txt')


if os.path.exists("./potential_loop/data.npz"):
    loaded_data = numpy.load("./potential_loop/data.npz")

    nsites_array = loaded_data["nsites_array"]
    nsources_array = loaded_data["nsources_array"]
    nmesh_array = loaded_data["nmesh_array"]
    noise_coeff_array = loaded_data["noise_coeff_array"]
    noise_level_array = loaded_data["noise_level_array"]
    finit_array = loaded_data["finit_array"]
    normgpinit_array = loaded_data["normgpinit_array"]
    ffinal_array = loaded_data["ffinal_array"]
    normgpfinal_array  = loaded_data["normgpfinal_array"]
    erroropt_array  = loaded_data["erroropt_array"]
    errorinit_array  = loaded_data["errorinit_array"]
    flagsol_array = loaded_data["flagsol_array"]
    iter_array = loaded_data["iter_array"]
    numevalf_array = loaded_data["numevalf_array"]
    CPU_time_array = loaded_data["CPU_time_array"]

    nsites_list = nsites_array.tolist()
    nsources_list = nsources_array.tolist()
    nmesh_list = nmesh_array.tolist()
    noise_coeff_list = noise_coeff_array.tolist()
    noise_level_list = noise_level_array.tolist()
    finit_list = finit_array.tolist()
    normgpinit_list = normgpinit_array.tolist()
    ffinal_list = ffinal_array.tolist()
    normgpfinal_list = normgpfinal_array.tolist()
    erroropt_list = erroropt_array.tolist()
    errorinit_list = errorinit_array.tolist()
    flagsol_list = flagsol_array.tolist()
    iter_list = iter_array.tolist()
    numevalf_list = numevalf_array.tolist()
    CPU_time_list = CPU_time_array.tolist()

    nsites_list.append(nsites)
    nsources_list.append(nsources)
    nmesh_list.append(nmesh)
    noise_coeff_list.append(noise_coeff)
    finit_list.append(finit)
    noise_level_list.append(noise_level)
    normgpinit_list.append(normgpinit)
    ffinal_list.append(ffinal)
    normgpfinal_list.append(normgpfinal)
    erroropt_list.append(erroropt)
    errorinit_list.append(errorinit)
    flagsol_list.append(flagsol)
    iter_list.append(iter)
    numevalf_list.append(numevalf)
    CPU_time_list.append(CPU_time)
else:
    nsites_list = [nsites]
    nsources_list = [nsources]
    nmesh_list = [nmesh]
    noise_coeff_list = [noise_coeff]
    finit_list = [finit]
    noise_level_list =[noise_level]
    normgpinit_list = [normgpinit]
    ffinal_list = [ffinal]
    normgpfinal_list = [normgpfinal]
    erroropt_list = [erroropt]
    errorinit_list = [errorinit]
    flagsol_list = [flagsol]
    iter_list = [iter]
    numevalf_list = [numevalf]
    CPU_time_list = [CPU_time]

nsites_array = numpy.array(nsites_list)
nsources_array = numpy.array(nsources_list)
nmesh_array = numpy.array(nmesh_list)
noise_coeff_array = numpy.array(noise_coeff_list)
finit_array = numpy.array(finit_list)
noise_level_array = numpy.array(noise_level_list)
normgpinit_array = numpy.array(normgpinit_list)
ffinal_array = numpy.array(ffinal_list)
normgpfinal_array = numpy.array(normgpfinal_list)
erroropt_array = numpy.array(erroropt_list)
errorinit_array = numpy.array(errorinit_list)
flagsol_array = numpy.array(flagsol_list)
iter_array = numpy.array(iter_list)
numevalf_array = numpy.array(numevalf_list)
CPU_time_array = numpy.array(CPU_time_list)

numpy.savez("./potential_loop/data.npz", nsites_array=nsites_array, nsources_array=nsources_array, nmesh_array=nmesh_array, noise_coeff_array=noise_coeff_array,
finit_array=finit_array,
noise_level_array=noise_level_array,
normgpinit_array=normgpinit_array,
ffinal_array=ffinal_array,
normgpfinal_array=normgpfinal_array, erroropt_array=erroropt_array,
errorinit_array=errorinit_array,
flagsol_array=flagsol_array, iter_array=iter_array, numevalf_array=numevalf_array, CPU_time_array=CPU_time_array)


print('End of main loop !!')
