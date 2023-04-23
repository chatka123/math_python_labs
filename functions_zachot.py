import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

def divided_difference(X,ind,por,f):
    X = X[ind-1:ind+por-1]
    Y = []
    for a in X:
        Y.append(f(a))
    value = []
    value2 = []
    for i in range(len(X)):
        difference = (Y[i+1] - Y[i])/(X[i+1]-X[i]) #сохраняем первый порядок по порядку
        value.append(difference)
        print(difference, ' - разделенная разность первого порядка')
        if i+2 == len(X):
            break

    order = 2
    if len(value) == 1:
        pass
        
    else:
        flag = True
        while flag:
            if flag == True:
                for i in range(len(X)):
                    difference = (value[1]-value[0])/(X[i+order]-X[i]) 
                    value2.append(difference)
                    print(difference, ' - разделенная разность порядка',order)
                    del value[0]
                    if len(value) == 1:
                        if len(value2) == 1:
                            flag = False
                            break
                        else:
                            del value[0]
                            order += 1
                            print('Следующий порядок = ', order)
                            break
            else:
                break


            if flag == True:
                for i in range(len(X)):
                    difference = (value2[1]-value2[0])/(X[i+order]-X[i])
                    value.append(difference)
                    print(difference, ' - разделенная разность порядка',order+1)
                    del value2[0]
                    if len(value2) == 1:
                        if len(value) ==1:
                            flag = False
                            break
                        else:
                            del value2[0]
                            order += 1
                            print('Следующий порядок = ', order+1)
                            break
            else:
                break

#_______________________________________________________________________________________________________________________________________________

def coeffs_newton(X,Y,f,X_funk,Y_funk):
    
    value = []
    value2 = []                 #пара 'рабочих' списка, в которые сохраняются разности разных порядков
    
    coeffs = [Y[0]] #коэффициенты будут сохраняться тут 
    count = 0   #счетчик будет следить за итерациями циклов чтобы брать как коэфф ньютона именно первую разность
    
    for i in range(len(X)):
        difference = (Y[i+1] - Y[i])/(X[i+1]-X[i]) #сохраняем первый порядок по порядку
        value.append(difference)
        print(difference, ' - разделенная разность первого порядка')
        if count == 0: 
            coeffs.append(difference)
            count = 0
        count += 1
        if i+2 == len(X):
            count = 0
            break

    order = 2  #счетчик который следит за порядком разности
    
    if len(value) == 1:
        pass
    else:
        flag = True
        while flag:
            if flag == True:
                for i in range(len(X)):
                    difference = (value[1]-value[0])/(X[i+order]-X[i])
                    value2.append(difference)
                    print(difference, ' - разделенная разность порядка',order)
                    del value[0]
                    if count == 0:
                        coeffs.append(difference)
                    count +=1 #увеличиваем счетчик чтобы следующие разности не попали в список коеффов
                    if len(value) == 1:
                        if len(value2) == 1:
                            flag = False
                            break
                        else:
                            del value[0]
                            order += 1
                            print('Следующий порядок = ', order)
                            count = 0 #обнуляем счетчик итераций чтобы далее работал корректно
                            break
            else:
                break


            if flag == True:
                for i in range(len(X)):
                    difference = (value2[1]-value2[0])/(X[i+order]-X[i])
                    value.append(difference)
                    print(difference, ' - разделенная разность порядка',order)
                    del value2[0]
                    if count == 0:
                        coeffs.append(difference)
                    count +=1 #увеличиваем счетчик чтобы следующие разности не попали в список коеффов
                    if len(value2) == 1:
                        if len(value) ==1:
                            flag = False
                            break
                        else:
                            del value2[0]
                            order += 1
                            print('Следующий порядок = ', order)
                            count = 0 #обнуляем счетчик итераций чтобы далее работал корректно
                            break
            else:
                break

    print('Coefficients of Newton`s polynomial ',coeffs)
    plot(X,Y,coeffs,f,X_funk,Y_funk)

#_______________________________________________________________________________________________________________________________________________
    
    
    
def plot(X,Y,coeffs,f,X_funk,Y_funk):
    print('Красный - функция')
    print('Синий - ньютоновский многочлен')
    plt.plot(X_funk,Y_funk,'r')
    Y_newton = y_newton(X,Y,coeffs)
    plt.plot(X,Y_newton,'b')
    plt.grid()
    plt.show()
    i = int(input('Относительно какой точки посчитать приращение функции и дифференциал? (индекс)'))
    if i <= len(X):
        x = sp.symbols('x')
        x0 = X[i-1]
        _f = f.subs(x,x0+0.1)-f.subs(x,x0)
        df = sp.diff(f,x,1).subs(x,x0)*0.1
        print('Индексу соответствует x = ',x0)    
        print('приращение функции = ', round(_f,4))
        print('дифференциал функции = ', round(df,4))
        print('их отношение = ', round(df/_f,4))
        tangent_line(f,X,X_funk,Y_funk)
    else:
        print('Точка с таким индексом не была найдена')
        

#_______________________________________________________________________________________________________________________________________________
        
def y_newton(X,Y,coeffs):
    y_var = coeffs[0]
    y_values = [coeffs[0]]
    for i in range(len(X)-1):
        x0 = X[i+1]
        for j in range(i+1):
            y_newton = coeffs[j+1]
            for k in range(j+1):
                y_newton *= (x0-X[k])
            y_var += y_newton
        y_values.append(y_var)
        y_var = coeffs[0]
    return y_values      
        
        

        
        
        

#___________________________________________________________________________________________________________________-
        
        
        
def y_newton_lose(X,Y,coeffs):
    '''Составляет уравнение многочлена ньютона'''
    x = np.linspace(X[0],X[-1],len(X))
    y_newton = coeffs[0]
    newton_coeffs = []
    for i in range(1,len(coeffs)):
        y_var = coeffs[i]
        print('Current coefficient: ', y_var)
        for j in range(i):
            y_var *= (x - X[j])
            newton_coeffs.append(y_var)
            
    y_newton += sum(newton_coeffs)
    return y_newton

#_______________________________________________________________________________________________________________________________________________



def tangent_line(f,xn_,X_funk,Y_funk):
    '''Строит нормаль к графику функции в каждой получаемой точке'''
    for i in range(len(xn_)):
        plt.plot(X_funk,Y_funk,'r')
        x = sp.symbols('x')
        f_dx = sp.diff(f,x,1)
        f0 = f.subs(x,xn_[i])
        f_dx = f_dx.subs(x,xn_[i])
        x = np.linspace(-Y_funk,Y_funk,len(X_funk))
        concern = f0+(1/f_dx)*(x-xn_[i])
        plt.plot(x,concern,'g')
        plt.grid()
        plt.show()

    
    
    
#_______________________________________________________________________________________________________________________________________________