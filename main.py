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

#процедура составляющая список уравнений для сплайна по массивам точек X,Y points(x[i],y[i])
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


def spline(X,Y):
    x=sympy.Symbol('x')
    a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3=sympy.symbols('a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3')

    nscA=a0+a1*x+a2*x**2+a3*x**3
    nscB=b0+b1*x+b2*x**2+b3*x**3
    nscC=c0+c1*x+c2*x**2+c3*x**3
    nscD=d0+d1*x+d2*x**2+d3*x**3

    nsc = [nscA, nscB, nscC, nscD]
    #X = [x1, x2, x3, x4, x5]
    #Y = [y1, y2, y3, y4, y5]
    eq = list_of_equations(X, Y, nsc)
    aa = [a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3]

    sol2 = sympy.solve(eq, aa)
    for i in range(len(aa)):
        aa[i] = sol2[aa[i]]

    #запишем коэффициенты в файл
    original_stdout = sys.stdout
    Coefficients = open('Coefficients', 'w')
    sys.stdout = Coefficients
    for i in range(16):
        if (i//4==0):
            letter='a'
        if (i//4==1):
            letter='b'
        if i//4==2:
            letter='c'
        if i//4==3:
            letter='d'
        print(letter+str(i%4)+'=')
        print(aa[i])

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
def list_of_symbols(n,j,letter):
    list_symb=[None]*n
    for i in range(n):
        stri=letter+str(j)+str(i)
        stri = sympy.Symbol(letter+str(j)+str(i))
        list_symb[i]=stri
    return list_symb



def main():
   x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5=sympy.symbols('x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5')

   x=sympy.Symbol('x')
   #a0,a1,a2,a3=polynom_with_boundary_conditions(x0,y0,x2,y2)
   #f=a0+a1*x+a2*x**2+a3*x**3

   #a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3=sympy.symbols(' a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3')

   X=[x1,x2,x3,x4,x5]
   Y=[y1,y2,y3,y4,y5]
   [a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3]=spline(X,Y)

   X_values=[0,1,2,3,4]
   Y_values=[0,1,-2,-1,-1]



   #нарисуем график
   aa=[a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3]
   for i in range(len(aa)):
       aa[i]=aa[i].subs({x1: X_values[0], y1: Y_values[0], x2: X_values[1], y2: Y_values[1], x3: X_values[2], y3: Y_values[2], x4: X_values[3],
                         y4: Y_values[3], x5: X_values[4], y5: Y_values[4]})

   #


   nsc = [[a0,a1,a2,a3], [b0,b1,b2,b3], [c0,c1,c2,c3], [d0,d1,d2,d3]]
   nsc = [[aa[0], aa[1], aa[2], aa[3]], [aa[4], aa[5], aa[6], aa[7]], [aa[8], aa[9], aa[10], aa[11]], [aa[12], aa[13], aa[14], aa[15]]]
   #X = [x11, x22, x33, x44, x55]
   draw_pieswize_f(nsc,X_values )



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
