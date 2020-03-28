import pylab as pl

def linear(x, y, dy):
    # LEAST-SQUARE FIT da Frodesen-Skjeggestad-Tofte
    # Assunzione: osservazioni indipendenti!
    Vy = pl.diag(pl.square(dy))
    Vy_inv = pl.inv(Vy)
    o = pl.ones(len(x))
    At = pl.asmatrix([x, o])
    A = pl.transpose(At)

    covar_m = pl.inv(pl.dot(At, pl.dot( Vy_inv, A)))
    V_invY = pl.dot(Vy_inv, y)
    AtVy_invY = pl.dot(At, V_invY)
    pars_v = pl.dot(covar_m, AtVy_invY.A1).A1
    dpars_v = pl.sqrt(pl.diagonal(covar_m))
    
    print('Fit lineare: minimi quadrati')  
    print('Parametri di fit (m,q): ',pars_v)
    print('Errori           (m,q): ',dpars_v)
    print('Matrice di covarianza:')
    print(covar_m)
    print('Correlazione     (m,q): ',covar_m[0,1]/(dpars_v[0]*dpars_v[1]))
    
    # Definisco la funzione
    def f(x,m,q):
        return (m*x + q)

    # Calcolo del Chi^2
    chi = ((y - f(x, pars_v[0], pars_v[1]))**2/dy**2).sum()
    ndof = len(x) - 2
    chi_red = chi/ndof
    
    chi_v = pl.array([chi, ndof])

    print('Gradi di libertà: ',ndof)
    print('Chiquadro    : ', chi)
    print('Chi^2 ridotto: ', chi_red) 

    return pars_v, dpars_v, covar_m, chi_v

def linear_qfixed(x, y, dy, q):
    # LEAST-SQUARE FIT da Frodesen-Skjeggestad-Tofte
    # Assunzione: osservazioni indipendenti!
    Vy = pl.diag(pl.square(dy))
    Vy_inv = pl.inv(Vy)
    o = pl.ones(len(x))
    At = pl.asmatrix([x])
    A = pl.transpose(At)

    covar_m = pl.inv(pl.dot(At, pl.dot( Vy_inv, A)))
    V_invY = pl.dot(Vy_inv, y)
    AtVy_invY = pl.dot(At, V_invY)
    pars_v = pl.dot(covar_m, AtVy_invY.A1).A1
    dpars_v = pl.sqrt(pl.diagonal(covar_m))
    
    print('Fit lineare: minimi quadrati')  
    print('Parametri di fit (m): ',pars_v)
    print('Errori           (m): ',dpars_v)
    print('Matrice di covarianza:')
    print(covar_m)
    #print('Correlazione     (m,q): ',covar_m[0,1]/(dpars_v[0]*dpars_v[1]))
    
    # Definisco la funzione
    def f(x,m,q):
        return (m*x + q)

    # Calcolo del Chi^2
    chi = ((y - f(x, pars_v[0], q))**2/dy**2).sum()
    ndof = len(x) - 1
    chi_red = chi/ndof
    
    chi_v = pl.array([chi, ndof])

    print('Gradi di libertà: ',ndof)
    print('Chiquadro    : ', chi)
    print('Chi^2 ridotto: ', chi_red) 

    return pars_v, dpars_v, covar_m, chi_v


def linear_mfixed(x, y, dy, m):
    # LEAST-SQUARE FIT da Frodesen-Skjeggestad-Tofte
    # Assunzione: osservazioni indipendenti!
    Vy = pl.diag(pl.square(dy))
    Vy_inv = pl.inv(Vy)
    o = pl.ones(len(x))
    At = pl.asmatrix([o])
    A = pl.transpose(At)

    covar_m = pl.inv(pl.dot(At, pl.dot( Vy_inv, A)))
    V_invY = pl.dot(Vy_inv, y)
    AtVy_invY = pl.dot(At, V_invY)
    pars_v = pl.dot(covar_m, AtVy_invY.A1).A1
    dpars_v = pl.sqrt(pl.diagonal(covar_m))
    
    print('Fit lineare: minimi quadrati')  
    print('Parametri di fit (q): ',pars_v)
    print('Errori           (q): ',dpars_v)
    print('Matrice di covarianza:')
    print(covar_m)
    #print('Correlazione     (m,q): ',covar_m[0,1]/(dpars_v[0]*dpars_v[1]))
    
    # Definisco la funzione
    def f(x,m,q):
        return (m*x + q)

    # Calcolo del Chi^2
    chi = ((y - f(x, m, pars_v[0]))**2/dy**2).sum()
    ndof = len(x) - 1
    chi_red = chi/ndof
    
    chi_v = pl.array([chi, ndof])

    print('Gradi di libertà: ',ndof)
    print('Chiquadro    : ', chi)
    print('Chi^2 ridotto: ', chi_red) 

    return pars_v, dpars_v, covar_m, chi_v
