import numpy as np
import sympy
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
import sys #для записи в файл
from sympy.plotting import (plot_parametric,plot_implicit,plot3d,plot3d_parametric_line,plot3d_parametric_surface)


#процедура вычисляющая коэффициента многочлена f(x)=a0+a1x+a2x**2+a3x**3
# для которого f(x0)=y0,f(x1)=y1,f'(x0)=0, f'(x1)=0
#эта процедура нигде не нужна
def polynom_with_boundary_conditions(x0,y0,x1,y1):
    d=x1-x0
    h=y0-y1
    a0=y0-3*h*x0**2/(d**2)-2*h*x0**3/(d**3)
    a1=6*h*x0/(d**2)+6*h*x0**2/(d**3)
    a2=-3*h/(d**2)-6*h*x0/(d**3)
    a3=2*h/(d**3)
    return a0,a1,a2,a3

#процедура составляющая список уравнений для 1-сплайна по массивам точек X,Y points(x[i],y[i])
def list_of_equations(X,Y,nsc):
    x = sympy.Symbol('x')
    l=len(X)
    n=l-1
    eq=[]
    eq0=nsc[0].subs(x,X[0])-Y[0]
    eq1=diff(nsc[0], x).subs(x, X[0])
    eq2=nsc[l-2].subs(x,X[l-1])-Y[l-1]
    eq3 = diff(nsc[l-2], x).subs(x, X[l-1])
    eq.append(eq0)
    eq.append(eq1)
    eq.append(eq2)
    eq.append(eq3)
    for i in range(1,n):
        eq1=nsc[i-1].subs(x, X[i]) - Y[i]
        eq2=nsc[i].subs(x, X[i]) - Y[i]
        eq3=diff(nsc[i-1], x).subs(x, X[i]) - diff(nsc[i], x).subs(x, X[i])
        eq4=diff(nsc[i-1], x, 2).subs(x, X[i]) - diff(nsc[i], x, 2).subs(x, X[i])
        eq.append(eq1)
        eq.append(eq2)
        eq.append(eq3)
        eq.append(eq4)
    return eq

#процудура, добавляющая к списку eq уравнения, что все коэффициенты многочлена f равны нулю
def polynom_eq_zerо_equations(f,eq,x):
    #x = sympy.Symbol('x')
    degf=degree(f,x)
    for i in range(degf+1):
        eqi=diff(f, x,i).subs(x, 0)
        eq.append(eqi)

#процедура составляющая список уравнений для 2-сплайна по массивам точек X,Y,Z
# Z=[[F(x0,y0),F(x1,y0),F(x2,y0)..F(xn,y0)],[],..[]]
def list_of_equations_for_bispline(X,Y,Z,nsc):
    x,y=sympy.symbols('x,y')
    lx=len(X)
    ly=len(Y)
    eq=[]
    #q=0
    #уравнения для значений в точках
    for i in range(lx-1):
        for j in range(ly-1):
            eq1=nsc[j][i].subs({x:X[i],y:Y[j]})-Z[j][i]
            eq2 = nsc[j][i].subs({x: X[i+1], y: Y[j]}) - Z[j][i+1]
            eq3 = nsc[j][i].subs({x: X[i], y: Y[j+1]}) - Z[j+1][i]
            eq4 = nsc[j][i].subs({x: X[i+1], y: Y[j+1]}) - Z[j+1][i+1]
            eq.append(eq1)
            eq.append(eq2)
            eq.append(eq3)
            eq.append(eq4)
            #q+=4
            #print(q,'ij',i,j)
    #уравнения для значений первых производных на границе прямоугольника
    for i in range(lx-1):
        f=diff(nsc[0][i],y).subs(y,Y[0])
        g = diff(nsc[ly - 2][i], y).subs(y, Y[ly - 1])
        #q+=degree(f,x)+degree(g,x)
        polynom_eq_zerо_equations(f,eq,x)
        polynom_eq_zerо_equations(g, eq,x)
        #print(q, 'i', i)
    for i in range(ly-1):
        f=diff(nsc[i][0],x).subs(x,X[0])
        g=diff(nsc[i][lx-2],x).subs(x,X[lx-1])
        polynom_eq_zerо_equations(f, eq, y)
        polynom_eq_zerо_equations(g, eq, y)
        #q += degree(f, x) + degree(g, x)
        #print(q, 'j', i)
    # уравнения для значений первых и вторых производных
    for i in range(1,lx-1):
        for j in range(ly-1):
            f=diff(nsc[j][i-1],x).subs(x,X[i])-diff(nsc[j][i],x).subs(x,X[i])
            g=diff(nsc[j][i-1],x,2).subs(x,X[i])-diff(nsc[j][i],x,2).subs(x,X[i])
            polynom_eq_zerо_equations(f, eq, y)
            polynom_eq_zerо_equations(g, eq, y)
    for j in range(1,ly-1):
        for i in range(lx-1):
            f=diff(nsc[j-1][i],y).subs(y,Y[j])-diff(nsc[j][i],y).subs(y,Y[j])
            g = diff(nsc[j - 1][i], y,2).subs(y, Y[j]) - diff(nsc[j][i], y,2).subs(y, Y[j])
            polynom_eq_zerо_equations(f, eq, x)
            polynom_eq_zerо_equations(g, eq, x)
    return eq

def spline(X,Y):
    n=len(X)
    if n!=len(Y):
        print('something goes wrong')
    x=sympy.Symbol('x')
    aa=list_of_2dsymbols(n-1,4,'a')
    aa_as_one_list=[]#нужен будет только как список переменных, относительно которых система уравнений
    for i in aa:
        aa_as_one_list.extend(i)

    nsc=[]
    for i in range(n-1):
        nsc.append(aa[i][0]+aa[i][1]*x+aa[i][2]*x**2+aa[i][3]*x**3)

    eq = list_of_equations(X, Y, nsc)
    sol2 = sympy.solve(eq, aa_as_one_list)
    for i in range(n-1):
        for j in range(4):
            aa[i][j]=sol2[aa[i][j]]

    #запишем коэффициенты в файл
    original_stdout = sys.stdout
    Coefficients = open('Coefficients', 'w')
    sys.stdout = Coefficients
    for i in range(n-1):
        for j in range(4):
            print('a'+str(i)+str(j)+'=')
            print(aa[i][j])
    sys.stdout = original_stdout
    Coefficients.close()
    return aa


# Z=[[F(x0,y0),F(x1,y0),F(x2,y0)..F(xn,y0)],[],..[]]
def bi_spline(X,Y,Z):
    x, y = sympy.symbols('x,y')
    n=len(X)
    m=len(Y)
    if m!=len(Z):
        print('something goes wrong')
    aa=list_of_4dsymbols(n-1,m-1,4,4,'a')
    aa_as_one_list=[]#нужен будет только как список переменных, относительно которых система уравнений
    for i in range(n-1):
        for j in range(m-1):
            for k in range(4):
                for l in range(4):
                    aa_as_one_list.append(aa[j][i][l][k])
    nsc=[] #создадим список многочленов для сплайна
    for i in range(m-1):
        nsci=[]
        for j in range(n-1):
            f=0
            for l in range(4):
                for k in range(4):
                    f=f+aa[i][j][k][l]*x**l*y**k
            nsci.append(f)
        nsc.append(nsci)

    eq=list_of_equations_for_bispline(X, Y, Z, nsc)
    sol = sympy.solve(eq, aa_as_one_list)

    for i in range(m-1):
        for j in range(n-1):
            for l in range(4):
                for k in range(4):
                    #aa[i][j][k][l]=1
                    aa[i][j][k][l]=sol[aa[i][j][k][l]]

    #запишем коэффициенты в файл
    original_stdout = sys.stdout
    BiCoefficients = open('BiCoefficients', 'w')
    sys.stdout = BiCoefficients
    for i in range(m - 1):
        for j in range(n - 1):
            for l in range(4):
                for k in range(4):
                    print('a' + str(i) + str(j) + str(k) + str(l) + '=')
                    print(aa[i][j][k][l])

    sys.stdout = original_stdout
    BiCoefficients.close()
    return aa


#nsc - список коэффициентов многочленов
#X - список точек, определяющий сегменты на которых определены функции nsc
def draw_pieswize_f(nsc,X):
    x = sympy.Symbol('x')
    l=len(nsc)
    p0=plot(nsc[0][0] + nsc[0][1] * x + nsc[0][2] * x ** 2 + nsc[0][3] * x ** 3,(x,X[0],X[1]),show=False)
    for i in range(1,l):
        p1=plot(nsc[i][0] + nsc[i][1] * x + nsc[i][2] * x ** 2 + nsc[i][3] * x ** 3,(x,X[i],X[i+1]),show=False)
        p0.extend(p1)
    p0.show()

#nsc - список функций
def draw_3d_pieswize_f(nsc,X,Y):
    x, y = sympy.symbols('x,y')
    m=len(nsc)+1
    n=len(nsc[0])+1
    p0=plot3d(nsc[0][0],(x,X[0],X[1]),
            (y,Y[0],Y[1]),show=False)
    for i in range(m-1):
        for j in range(n-1):
            p1=plot3d(nsc[i][j],(x,X[j],X[j+1]),
            (y,Y[i],Y[i+1]),show=False)
            p0.extend(p1)
    p0.show()

#процедура, создающая список символьных переменных
def list_of_symbols(n,letter):
    #letter - буква, используемая для переменной
    list_symb=[None]*n
    for i in range(n):
        stri = sympy.Symbol(letter+str(i))
        list_symb[i]=stri
    return list_symb

#процедуры, создающая список списков символьных переменных
def list_of_2dsymbols(n,m,letter):
    list_symb=[None]*n
    for i in range(n):
        stri = list_of_symbols(m,letter+str(i))
        list_symb[i]=stri
    return list_symb


def list_of_3dsymbols(n,m,k,letter):
    list_symb=[None]*m
    for i in range(m):
        stri = list_of_2dsymbols(n,k,letter+str(i))
        list_symb[i]=stri
    return list_symb

def list_of_4dsymbols(n,m,k,l,letter):
    list_symb=[None]*m
    for i in range(m):
        stri = list_of_3dsymbols(l,n,k,letter+str(i))
        list_symb[i]=stri
    return list_symb

#процудура создающая список многочленов по списку коэффициентов
def create_bi_nsc(BiSplineCoeff):
    x, y = sympy.symbols('x,y')
    m=len(BiSplineCoeff)+1
    n=len(BiSplineCoeff[0])+1
    print('m,n',m,n)
    nsc = []
    for j in range(m - 1):
        nscj = []
        for i in range(n - 1):
            f = 0
            for ii in range(4):
                for jj in range(4):
                    f += BiSplineCoeff[j][i][jj][ii] * x ** ii * y ** jj
            nscj.append(f)
        nsc.append(nscj)
    return nsc

def main():

   x, y = sympy.symbols('x,y')
   #a0,a1,a2,a3=polynom_with_boundary_conditions(x0,y0,x2,y2)
   #f=a0+a1*x+a2*x**2+a3*x**3

   n=4

   X=list_of_symbols(n,'x')
   Y=list_of_symbols(n,'y')

   SplineCoeff = spline(X, Y)

   X_values=[-1,1,2,3]
   Y_values=[0,1,-2,-1]

   #подставим значения в выражения для коэффициентов
   sc_as_one_list = []
   for i in SplineCoeff:
       sc_as_one_list.extend(i)


   for i in range(len(sc_as_one_list)):
       for j in range(n):
           sc_as_one_list[i]=sc_as_one_list[i].subs({X[j]:X_values[j],Y[j]:Y_values[j]})

   # нарисуем график
   nsc=[]
   for i in range(n-1):
       nsc.append([sc_as_one_list[4*i+0],sc_as_one_list[4*i+1],sc_as_one_list[4*i+2],sc_as_one_list[4*i+3]])

   draw_pieswize_f(nsc,X_values )

   #БИСПЛАЙН
   n=3
   m=3
   X = list_of_symbols(n, 'x')
   Y = list_of_symbols(m, 'y')
   Z=list_of_2dsymbols(m,n,'z')


   BiSplineCoeff = bi_spline(X, Y,Z)
   print('BiSplineCoeff',BiSplineCoeff)

   X_values = [0,1,3]
   Y_values = [0, 1,2]
   Z_values=[[0,0,0],[1,3,2],[1,2,4]]
   # подставим значения в выражения для коэффициентов
   bsc_as_one_list = []
   for i in range(n - 1):
       for j in range(m - 1):
           for k in range(4):
               for l in range(4):
                   for ii in range(n):
                       for jj in range(m):
                           BiSplineCoeff[j][i][k][l]=BiSplineCoeff[j][i][k][l].subs(
                               {X[ii]: X_values[ii], Y[jj]: Y_values[jj], Z[jj][ii]: Z_values[jj][ii]})
   print(BiSplineCoeff)

   # нарисуем график
   nsc = create_bi_nsc(BiSplineCoeff)
   print('nsc',nsc)

   draw_3d_pieswize_f(nsc, X_values, Y_values)

   #проверка
   #разница нужных и полученных значений в точках
   for i in range(n-1):
       for j in range(m-1):
           d1=nsc[j][i].subs({x:X_values[i],y:Y_values[j]})-Z_values[j][i]
           d2 = nsc[j][i].subs({x:X_values[i+1],y:Y_values[j]})-Z_values[j][i+1]
           d3 = nsc[j][i].subs({x:X_values[i],y:Y_values[j+1]})-Z_values[j+1][i]
           d4 = nsc[j][i].subs({x:X_values[i+1],y:Y_values[j+1]})-Z_values[j+1][i+1]
           print('di',i,j,d1,d2,d3,d4)

       # разница нужных и полученных значений производных
   for i in range(1,n - 1):
       for j in range(m - 1):
           d1 = diff(nsc[j][i-1],x).subs({x: X_values[i]}) - diff(nsc[j][i],x).subs({x: X_values[i]})
           print('di', i, j, d1)
   for j in range(1,m - 1):
       for i in range(n - 1):
           d1 = diff(nsc[j-1][i],y).subs({y: Y_values[j]}) - diff(nsc[j][i],y).subs({y: Y_values[j]})
           print('di', i, j, d1)





# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
