import numpy as np
import sympy
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
import sys #для записи в файл


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

#процедура составляющая список уравнений для 2-сплайна по массивам точек X,Y,Z
# Z=[[F(x0,y0),F(x0,y1),F(x0,y2)..F(x0,yn)],[],..[]]
def list_of_equations_for_bispline(X,Y,Z,nsc):
    x,y=sympy.symbols('x,y')
    lx=len(X)
    ly=len(Y)
    eq=[]
    #уравнения для значений в точках
    for i in range(lx):
        for j in range(ly):
            eqij=nsc[j][i].subs({x:X[i],y:Y[j]})-Z[i][j]
            eq.append(eqij)
    #уравнения для значений первых производных на границе прямоугольника
    for i in range(lx-1):
        eq1=diff(nsc[0][i],x).subs(x,X[0])
        eq2=diff(nsc[ly-2][i],x).subs(x,X[lx-1])
        eq.append(eq1)
        eq.append(eq2)
    for i in range(ly-1):
        eq1=diff(nsc[i][0],y).subs(y,Y[0])
        eq2=diff(nsc[i][lx-2],y).subs(y,Y[ly-1])
        eq.append(eq1)
        eq.append(eq2)
    # уравнения для значений первых и вторых производных
    for i in range(1,lx-1):
        for j in range(ly-1):
            eq1=diff(nsc[j][i-1],x).subs(x,X[i])-diff(nsc[j][i],x).subs(x,X[i])
            eq2=diff(nsc[j][i-1],x,2).subs(x,X[i])-diff(nsc[j][i],x,2).subs(x,X[i])
            eq.append(eq1)
            eq.append(eq2)
    for j in range(1,ly-1):
        for i in range(lx-1):
            eq1=diff(nsc[j-1][i],y).subs(y,Y[j])-diff(nsc[j][i],y).subs(y,Y[j])
            eq2 = diff(nsc[j - 1][i], y,2).subs(y, Y[j]) - diff(nsc[j][i], y,2).subs(y, Y[j])
            eq.append(eq1)
            eq.append(eq2)
    return eq

def spline(n,X,Y):
    x=sympy.Symbol('x')
    #a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3=sympy.symbols('a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3')#replase
    aa=list_of_2dsymbols(n-1,4,'a')#new
    aa_as_one_list=[]
    for i in aa:
        aa_as_one_list.extend(i)

   # nscA=a0+a1*x+a2*x**2+a3*x**3#replase
    #nscB=b0+b1*x+b2*x**2+b3*x**3#replase
    #nscC=c0+c1*x+c2*x**2+c3*x**3#replase
    #nscD=d0+d1*x+d2*x**2+d3*x**3#replace
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

#nsc - список коэффициентов многочленов
def draw_pieswize_f(nsc,X):
    x = sympy.Symbol('x')
    l=len(nsc)
    p0=plot(nsc[0][0] + nsc[0][1] * x + nsc[0][2] * x ** 2 + nsc[0][3] * x ** 3,(x,X[0],X[1]),show=False)
    for i in range(1,l):
        p1=plot(nsc[i][0] + nsc[i][1] * x + nsc[i][2] * x ** 2 + nsc[i][3] * x ** 3,(x,X[i],X[i+1]),show=False)
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

#процедура, создающая список списков символьных переменных
def list_of_2dsymbols(n,m,letter):
    list_symb=[None]*n
    for i in range(n):
        stri = list_of_symbols(m,letter+str(i))
        list_symb[i]=stri
    return list_symb

def main():

   x=sympy.Symbol('x')
   #a0,a1,a2,a3=polynom_with_boundary_conditions(x0,y0,x2,y2)
   #f=a0+a1*x+a2*x**2+a3*x**3

   n=4

   X=list_of_symbols(n,'x')
   Y=list_of_symbols(n,'y')

   SplineCoeff = spline(n, X, Y)

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
   #Z=[[z00,z10,z20,z30,z40,z50]]



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
