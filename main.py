import numpy as np
import sympy
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
import sys #для записи в файл


def spline(x_ax,y_ax):
    #x,y - массивы точек y[i] - значение в точке x[i]
    #x=sympy.Symbol('x')
    #x=sympy.symbols('x')
    x=Symbol('x')
    if len(x_ax)!=len(y_ax):
        print('length wrong')
    else:
        g=x+1+sympy.sin(x)
        f=sympy.zeros(len(x_ax))
        for i in range(len(x_ax)-1):
            f[i]=1
            #f[i]=(y_ax[i+1]-y_ax[i])/(x_ax[i+1]-x_ax[i])*x+(x_ax[i+1]*y_ax[i]-x_ax[i]*y_ax[i+1])/(x_ax[i+1]-x_ax[i])
    return f,g

#процедура вычисляющая коэффициента многочлена f(x)=a0+a1x+a2x**2+a3x**3
# для которого f(x0)=y0,f(x1)=y1,f'(x0)=0, f'(x1)=0
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


def spline(X,Y,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5):
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

    #a0 = sol2[a0]

    #запишем коэффициенты в файл
    original_stdout = sys.stdout
    FileGraphMatrix = open('Coefficients', 'w')
    sys.stdout = FileGraphMatrix
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
    FileGraphMatrix.close()

    return aa


def long_subs(a,X,Y,X_values,Y_values):
    for i in range(len(X)):
        a.subs((X[i]), X_values[i])
        a.subs(str(Y[i]),Y_values[i])
    aa=a
    return aa

#nsc - список коэффициентов многочленов
def draw_pieswize_f(nsc,X):
    x = sympy.Symbol('x')
    l=len(nsc)
    p=[None]*l
    p[0]=plot(nsc[0][0] + nsc[0][1] * x + nsc[0][2] * x ** 2 + nsc[0][3] * x ** 3,(x,X[0],X[1]),show=False)
    for i in range(1,l):
        p[1]=plot(nsc[i][0] + nsc[i][1] * x + nsc[i][2] * x ** 2 + nsc[i][3] * x ** 3,(x,X[i],X[i+1]),show=False)
        p[0].extend(p[1])

    p[0].show()

#процедура, создающая список символьных переменных
def list_of_symbols(n,j,letter):
    list_symb=[None]*n
    for i in range(n):
        stri=letter+str(j)+str(i)
        stri = sympy.Symbol(letter+str(j)+str(i))
        list_symb[i]=stri
    return list_symb

#процедура, создающая список символьных переменных
def list_of_symbols4(n,letter):
    list_symb4=[None]*n
    for i in range(n):
        list_symb4[i]=list_of_symbols(4,i,letter)
    return list_symb4

def main():
   #print('hello')
   #x=np.array([0,1,2])
   #y=np.array([0,1,2])
   #f=spline(x,y)
   #print(f)

   x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5=sympy.symbols('x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5')
   x=sympy.Symbol('x')
   a0,a1,a2,a3=polynom_with_boundary_conditions(x0,y0,x2,y2)
   print('koefficienti',a0,a1,a2,a3)
   f=a0+a1*x+a2*x**2+a3*x**3
   Df=a1+2*a2*x+3*a3*x**2
   print('f',f)
   print('fx0',sympy.simplify(f.subs(x,x0)))
   print('fx2', sympy.simplify(f.subs(x, x2)))
   print('dfx0', sympy.simplify(Df.subs(x, x0)))
   print('dfx2', sympy.simplify(Df.subs(x, x2)))


   a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3=sympy.symbols(' a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3')
   n=5
   X=list_of_symbols(n,0,'x')
   Y=list_of_symbols(n,0,'y')
   print('X',X)
   print('Y',Y)
   X=[x1,x2,x3,x4,x5]
   Y=[y1,y2,y3,y4,y5]
   [a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3]=spline(X,Y,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5)

   #a0, a1, a2, a3 = polynom_with_boundary_conditions(x1, y1, x2, y2)
   #b0, b1, b2, b3 = polynom_with_boundary_conditions(x2, y2, x3, y3)
   #c0, c1, c2, c3 = polynom_with_boundary_conditions(x3, y3, x4, y4)
   #d0, d1, d2, d3 = polynom_with_boundary_conditions(x4, y4, x5, y5)


   x11=0
   x22=1
   x33=2
   x44=3
   x55=4
   y11=0
   y22=1
   y33=-2
   y44=3
   y55=-1
   X_values=[0,1,2,3,4]
   Y_values=[0,1,-2,3,-1]
   #a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3=long_subs(a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3,
   #           x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,
   #           x11,y11,x22,y22,x33,y33,x44,y44,x55,y55)

   nscA = a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3
   nscB = b0 + b1 * x + b2 * x ** 2 + b3 * x ** 3
   nscC = c0 + c1 * x + c2 * x ** 2 + c3 * x ** 3
   nscD = d0 + d1 * x + d2 * x ** 2 + d3 * x ** 3

   #nsc = [nscA, nscB, nscC, nscD]
   #X = [x1, x2, x3, x4, x5]
   #Y = [y1, y2, y3, y4, y5]
   #eq = list_of_equations(X, Y, nsc)
  # print('eq', eq)
   #print('lenEq', len(eq))
   #a0,a1,a2,a3=polynom_with_boundary_conditions(x0,y0,x2,y2)#delete

   #нарисуем график

   a0 = a0.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   print('a0', a0)
   a0=long_subs(a0,X,Y,X_values,Y_values)

   print('a0',a0)
   a1 = long_subs(a1, X, Y, X_values, Y_values)
   a2 = long_subs(a2, X, Y, X_values, Y_values)
   a3 = long_subs(a3, X, Y, X_values, Y_values)
   b0 = long_subs(b0, X, Y, X_values, Y_values)
   b1 = long_subs(b1, X, Y, X_values, Y_values)
   b2 = long_subs(b2, X, Y, X_values, Y_values)
   b3 = long_subs(b3, X, Y, X_values, Y_values)
   c0=long_subs(c0,X,Y,X_values,Y_values)
   c1 = long_subs(c1, X, Y, X_values, Y_values)
   c2 = long_subs(c2, X, Y, X_values, Y_values)
   c3 = long_subs(c3, X, Y, X_values, Y_values)
   d0 = long_subs(d0, X, Y, X_values, Y_values)
   d1 = long_subs(d1, X, Y, X_values, Y_values)
   d2 = long_subs(d2, X, Y, X_values, Y_values)
   d3 = long_subs(d3, X, Y, X_values, Y_values)


   aa=[a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3]
   for i in aa:
       i=i.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   a1 = a1.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   a2 = a2.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   a3 = a3.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   b0 = b0.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   b1 = b1.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   b2 = b2.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   b3 = b3.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   c0 = c0.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   c1 = c1.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   c2 = c2.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   c3 = c3.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})

   d0 = d0.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   d1 = d1.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   d2 = d2.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})
   d3 = d3.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})

   #nscA = a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3
   print('a0',a0)
   ##plot(a0 + a1 * x+a2*x**2+a3*x**3,b0 + b1 * x+b2*x**2+b3*x**3,
       # c0 + c1 * x+c2*x**2+c3*x**3,d0 + d1 * x+d2*x**2+d3*x**3,(x,x11,x55),ylim=(-2,5))
   ##plot(a0+1 + a1 * x + a2 * x ** 2, (x, 0, 1))

   nsc = [[a0,a1,a2,a3], [b0,b1,b2,b3], [c0,c1,c2,c3], [d0,d1,d2,d3]]
   X = [x11, x22, x33, x44, x55]
   draw_pieswize_f(nsc,X )
   #p1=   plot(a0 + a1 * x+a2*x**2+a3*x**3,(x,x11,x22),show=False)
  # p1.line_color='cyan'
  # p2 = plot(b0 + b1 * x + b2 * x ** 2 + b3 * x ** 3, (x, x22, x33), show=False)
  # p2.line_color = 'red'
  # p3 = plot(c0 + c1 * x + c2 * x ** 2 + c3 * x ** 3, (x, x33, x44), show=False,line_color='red')
  # p4 = plot(d0 + d1 * x + d2 * x ** 2 + d3 * x ** 3, (x, x44, x55), show=False,line_color='blue')

   #p1.extend(p2)
  # p1.extend(p3)
  # p1.extend(p4)
  # p1.show()

   list_o=list_of_symbols(5,0,'x')
   print(list_o)
   list1=list_of_symbols4(5,'a')
   print(list1)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
