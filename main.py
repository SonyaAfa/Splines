import numpy as np
import sympy
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import scatter
import sys #для записи в файл

#natural cubic spiline
def ncs(Nodes):
    n=len(Nodes)-1
    a=np.zeros(n+1)
    j=0
    for i in Nodes:
        a[j]=i[1]
        j+=1
    b = np.zeros(n)
    d=np.zeros(n)
    h=np.zeros(n)
    for i in range(n):
        h[i]=Nodes[i+1][0]-Nodes[i][0]
    alpha=np.zeros(n)
    for i in range(1,n,1):
        alpha[i]=3/h[i]*(a[i+1]-a[i])-3/h[i-1]*(a[i]-a[i-1])
    c=np.zeros(n+1)
    l = np.zeros(n + 1)
    mu = np.zeros(n + 1)
    z = np.zeros(n + 1)
    l[0]=1
    mu[0]=0
    z[0]=0
    for i in range(1, n, 1):
        l[i]=2*(Nodes[i+1][0]-Nodes[i-1][0])-h[i-1]*mu[i-1]
        mu[i]=h[i]/l[i]
        z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i]
    l[n]=1
    z[n]=0
    c[n]=0
    for j in range(n-1,-1,-1):
        c[j]=z[j]-mu[j]*c[j+1]
        b[j]=(a[j+1]-a[j])/h[j]-(h[j]*(c[j+1]+2*c[j]))/3
        d[j]=(c[j+1]-c[j])/(3*h[j])
    output_set=np.zeros((n,5))
    for i in range(n):
        output_set[i][0]=a[i]
        output_set[i][1]=b[i]
        output_set[i][2] = c[i]
        output_set[i][3] = d[i]
        output_set[i][4] = Nodes[i][0]
    return(output_set)

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


#
def spline(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5):
    x=sympy.Symbol('x')
    a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3=sympy.symbols('a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3')

    nscA=a0+a1*x+a2*x**2+a3*x**3
    nscB=b0+b1*x+b2*x**2+b3*x**3
    nscC=c0+c1*x+c2*x**2+c3*x**3
    nscD=d0+d1*x+d2*x**2+d3*x**3

    #sol=sympy.solve([a0*x1-1,b0*x2-3],[a0,b0])
    sol2=sympy.solve([nscA.subs(x,x1)-y1,
                      diff(nscA,x).subs(x,x1),
                      nscA.subs(x,x2)-y2,
                      nscB.subs(x, x2) - y2,
                      diff(nscA,x).subs(x,x2)-diff(nscB,x).subs(x,x2),
                      diff(nscA,x,2).subs(x,x2)-diff(nscB,x,2).subs(x,x2),

                      nscB.subs(x, x3) - y3,
                      nscC.subs(x, x3) - y3,
                      diff(nscB, x).subs(x, x3) - diff(nscC, x).subs(x, x3),
                      diff(nscB, x, 2).subs(x, x3) - diff(nscC, x, 2).subs(x, x3),

                      nscC.subs(x, x4) - y4,
                      nscD.subs(x, x4) - y4,
                      diff(nscC, x).subs(x, x4) - diff(nscD, x).subs(x, x4),
                      diff(nscC, x, 2).subs(x, x4) - diff(nscD, x, 2).subs(x, x4),

                      nscD.subs(x, x5) - y5,
                      diff(nscD, x).subs(x, x5)
                      ],[a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3])
    #print('so',sol)
    print('sol2',sol2)
    print('sol2a0', sol2[a0])
    a0=sol2[a0]
    a1 = sol2[a1]
    a2 = sol2[a2]
    a3 = sol2[a3]

    b0 = sol2[b0]
    b1 = sol2[b1]
    b2 = sol2[b2]
    b3 = sol2[b3]

    c0 = sol2[c0]
    c1 = sol2[c1]
    c2 = sol2[c2]
    c3 = sol2[c3]

    d0 = sol2[d0]
    d1 = sol2[d1]
    d2 = sol2[d2]
    d3 = sol2[d3]
    original_stdout = sys.stdout
    FileGraphMatrix = open('Coefficients', 'w')
    sys.stdout = FileGraphMatrix

    print('a0=')
    print(a0)#+ '\n')
    print('a1=')
    print(a1)#+ '\n') + '\n')
    print('a2=')
    print(a2)# + '\n')
    print('a3=')
    print( a3 )#+ '\n')
    print('b0=')
    print( b0 )#+ '\n')
    print('b1=')
    print( b1 )#+ '\n')
    print('b2=')
    print( b2 )#+ '\n')
    print('b3=')
    print( b3 )#+ '\n')
    print('c0=')
    print( c0 )#+ '\n')
    print('c1=')
    print( c1 )#+ '\n')
    print('c2=')
    print( c2 )#+ '\n')
    print('c3=')
    print( c3 )#+ '\n')
    print('d0=')
    print( d0 )#+ '\n')
    print('d1=')
    print( d1 )#+ '\n')
    print('d2=')
    print( d2 )#+ '\n')
    print('d3=')
    print( d3 )#+ '\n')

    sys.stdout = original_stdout
    FileGraphMatrix.close()

    return a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3

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


def spline_new(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5):
    x=sympy.Symbol('x')
    a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3=sympy.symbols('a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3')

    nscA=a0+a1*x+a2*x**2+a3*x**3
    nscB=b0+b1*x+b2*x**2+b3*x**3
    nscC=c0+c1*x+c2*x**2+c3*x**3
    nscD=d0+d1*x+d2*x**2+d3*x**3

    nsc = [nscA, nscB, nscC, nscD]
    X = [x1, x2, x3, x4, x5]
    Y = [y1, y2, y3, y4, y5]
    eq = list_of_equations(X, Y, nsc)

    sol2=sympy.solve(eq,[a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3])
    #print('so',sol)
    print('sol2',sol2)
    print('sol2a0', sol2[a0])
    a0=sol2[a0]
    a1 = sol2[a1]
    a2 = sol2[a2]
    a3 = sol2[a3]

    b0 = sol2[b0]
    b1 = sol2[b1]
    b2 = sol2[b2]
    b3 = sol2[b3]

    c0 = sol2[c0]
    c1 = sol2[c1]
    c2 = sol2[c2]
    c3 = sol2[c3]

    d0 = sol2[d0]
    d1 = sol2[d1]
    d2 = sol2[d2]
    d3 = sol2[d3]
    original_stdout = sys.stdout
    FileGraphMatrix = open('Coefficients', 'w')
    sys.stdout = FileGraphMatrix

    print('a0=')
    print(a0)#+ '\n')
    print('a1=')
    print(a1)#+ '\n') + '\n')
    print('a2=')
    print(a2)# + '\n')
    print('a3=')
    print( a3 )#+ '\n')
    print('b0=')
    print( b0 )#+ '\n')
    print('b1=')
    print( b1 )#+ '\n')
    print('b2=')
    print( b2 )#+ '\n')
    print('b3=')
    print( b3 )#+ '\n')
    print('c0=')
    print( c0 )#+ '\n')
    print('c1=')
    print( c1 )#+ '\n')
    print('c2=')
    print( c2 )#+ '\n')
    print('c3=')
    print( c3 )#+ '\n')
    print('d0=')
    print( d0 )#+ '\n')
    print('d1=')
    print( d1 )#+ '\n')
    print('d2=')
    print( d2 )#+ '\n')
    print('d3=')
    print( d3 )#+ '\n')

    sys.stdout = original_stdout
    FileGraphMatrix.close()

    return a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3


def long_subs(a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3,
              x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,
              x11,y11,x22,y22,x33,y33,x44,y44,x55,y55):
    list=[a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3]
    for i in list:
        i = i.subs({x1: x11, y1: y11, x2: x22, y2: y22, x3: x33, y3: y33, x4: x44, y4: y44, x5: x55, y5: y55})

    return a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3




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
   a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3=spline_new(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5)

   #a0, a1, a2, a3 = polynom_with_boundary_conditions(x1, y1, x2, y2)
   #b0, b1, b2, b3 = polynom_with_boundary_conditions(x2, y2, x3, y3)
   #c0, c1, c2, c3 = polynom_with_boundary_conditions(x3, y3, x4, y4)
   #d0, d1, d2, d3 = polynom_with_boundary_conditions(x4, y4, x5, y5)

   #f = a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3
   #Df = a1 + 2 * a2 * x + 3 * a3 * x ** 2
   #print('f', f)
   #print('fx1', sympy.simplify(f.subs(x, x1)))
   #print('fx2', sympy.simplify(f.subs(x, x2)))
   #print('dfx1', sympy.simplify(Df.subs(x, x1)))
   #print('dfx2', sympy.simplify(Df.subs(x, x2)))
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
   p1=   plot(a0 + a1 * x+a2*x**2+a3*x**3,(x,x11,x22),show=False)
   p1.line_color='cyan'
   p2 = plot(b0 + b1 * x + b2 * x ** 2 + b3 * x ** 3, (x, x22, x33), show=False)
   p2.line_color = 'red'
   p3 = plot(c0 + c1 * x + c2 * x ** 2 + c3 * x ** 3, (x, x33, x44), show=False,line_color='red')
   p4 = plot(d0 + d1 * x + d2 * x ** 2 + d3 * x ** 3, (x, x44, x55), show=False,line_color='blue')

   p1.extend(p2)
   p1.extend(p3)
   p1.extend(p4)
   p1.show()



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
