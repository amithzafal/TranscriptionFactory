#!/usr/bin/env python
# coding: utf-8

# In[27]:


import lumicks.pylake as lk
import matplotlib.pyplot as plt
import os
from shapely.geometry import Polygon
import math
import numpy as np
from scipy.special import erf
from scipy.optimize import curve_fit
import scipy.special
import scipy.optimize
import numpy as np
import colour as co
try:
    os.chdir("D:\Project DNA")
    print(1)
except:
    os.chdir("F:\Project DNA")
    print(2)

def WLC_function(x,Lc,Lp,f_offset):
    
    #return(Lc*(1-1/2*(4.11/(F+offset)/Lp)**(1/2)+(F+offset)/Sm))
    return(4.11/Lp*(1/4*(1-x/Lc)**(-2)-1/4+x/Lc)+f_offset)
f = 78125

def tri(A,B,C):
    for i in range(len(A)):
        for j in range(len(A)-1):
            if A[j+1] < A[j]:
                A[j+1], A[j] = A[j], A[j+1]
                B[j+1], B[j] = B[j], B[j+1]
                C[j+1], C[j] = C[j], C[j+1]
    return(A,B,C)

def fiit(A,B):
    from sklearn.linear_model import LinearRegression
    x = np.array(A)
    y = np.array(B)
    x = x.reshape(-1,1)
    y = y.reshape(-1,1)
    reg = LinearRegression(normalize=True)
    reg.fit(x,y)
    a = reg.coef_
    b = reg.intercept_
    absi = np.linspace(min(A),max(A),1000)
    return(absi,a,b)
    

def traitement(distance, force):
    d_temp = []
    f_temp = []
    for i in range(len(distance)):
        if distance[i]>1:
            d_temp.append(distance[i])
            f_temp.append(force[i])
    list_d = []
    list_f = []
    mini = min(d_temp)
    maxi = max(d_temp)
    a = 0
    b = 0
    d = []
    f = []
    for i in range(len(d_temp)):
        if (d_temp[i]-mini)/(maxi-mini)>0.9:
            if a == 0 and b == 0 :
                pass
            elif a == 1 and b == 1 :
                a = 0
                d.append(d_temp[i])
                f.append(f_temp[i])
            else:
                d.append(d_temp[i])
                f.append(f_temp[i])
        elif (d_temp[i]-mini)/(maxi-mini)<0.1:
            if a == 0 and b == 1 :
                b = 0
                a = 1
                list_d.append(d)
                list_f.append(f)
                #plt.plot(d,f)
                #plt.title("traitement")
                #plt.show()
                d = []
                f = []
            if a == 0 and b == 0 : 
                a = 1
            if a == 1 and b == 0 :
                pass
        else:
            if a == 0 and b == 0 :
                pass
            if a == 1 and b == 0 :
                b = 1
                d.append(d_temp[i])
                f.append(f_temp[i])
            if a == 1 and b == 1 : 
                d.append(d_temp[i])
                f.append(f_temp[i])
            if a == 0 and b == 1:
                d.append(d_temp[i])
                f.append(f_temp[i])
    return(list_d,list_f)

def erfunc(x, mFL, a, b):
    return mFL*erf((x-a)/(b*np.sqrt(2)))


color = ["blue",'green',"orange",'red']

def cut(distance,force): #CUT POUR EVITER L'HYSTERESIS QUAND ON WLC
    distance_2 = []
    force_2 = []    
    i = 0
    f = 0.8*max(force)
    while force[i]<f:
        distance_2.append(distance[i])
        force_2.append(force[i])
        i+=1
    return(distance_2,force_2)
   
    
    
    
def fit_WLC(force,distance,path):
    d,f = cut(distance,force)
    plt.plot(d,f)
    plt.show()
    #if input('Save [y,n] : ') == "y":
    #    file = open("Data_lp.txt",'a')
    #    for i in range(len(d)):
    #        file.writelines(f'{d[i]};{f[i]}\n')
    #    file.close()

    min_f = min(f)
    j = f.index(min_f)

    force_cut = f
    distance_cut = d
   
    for i in range(len(f)):
        if i<j:
            force_cut.append(f[j]-(f[i]-f[j]))
            distance_cut.append(d[i])
        else:
            force_cut.append(f[i])
            distance_cut.append(d[i])

    f = []
    d = []
    for i in range(len(force_cut)):
        if force_cut[i]<10 and distance_cut[i]>10:
            f.append(force_cut[i])
            d.append(distance_cut[i])
    debut = d[0]
    distance = d#[dist-debut+10.5 for dist in d]
    
    params, extras = scipy.optimize.curve_fit(WLC_function, np.array(distance),np.array(f),p0=[16.5,56,0]) # Marker 1
    fig  = plt.figure(1,figsize = (16, 12))
    plt.plot(distance,f)
    plt.plot(distance,WLC_function(np.array(f),*params))
    plt.show()
    print(f"Lc = {float(params[0])} f_offset = {float(params[2])}, Lp = {float(params[1])}")
    return(float(params[1]),float(params[0]),float(params[2])) #
    
  
    
    #model = lk.inverted_odijk("DNA").subtract_independent_offset() + lk.force_offset("DNA")
    #fit = lk.FdFit(model)
    #fit.add_data(path, f, d)
    #fit['DNA/d_offset'].value = 0
    #fit['DNA/d_offset'].lower_bound = -1
    #fit['DNA/d_offset'].upper_bound = 1
    #fit['DNA/f_offset'].value = 0
    #fit['DNA/f_offset'].lower_bound = -1
    #fit['DNA/f_offset'].upper_bound = 1
    #fit['DNA/Lp'].value = 50
    #fit['DNA/Lp'].lower_bound = 0
    #fit['DNA/Lp'].upper_bound = 1000
    #fit['DNA/Lc'].value = 16
    #fit['DNA/Lc'].lower_bound = 0
    #fit['DNA/Lc'].upper_bound = 100
    #fit['DNA/St'].value = 1500
    #fit['DNA/St'].lower_bound = 0
    #fit['DNA/St'].upper_bound = 2000
    #fit.fit()
    #return(fit['DNA/Lp'].value,fit['DNA/Lc'].value,fit)
#fit['DNA/St'].value,fit['DNA/d_offset'].value,fit['DNA/f_offset'].value) 

#f(d) = argmin[f](norm(DNA.Lc * (1 - (1/2)*sqrt(kT/(f*DNA.Lp)) + f/DNA.St)-(d - DNA.d_offset))) + DNA.f_offset
  
def Moy_and_Var(values):
    moy_total = moy(values[1])
    abscisse = []
    list_moy = []
    temp = []
    list_var = []
    a = values[0][0]
    print(values)
    for i in range(len(values[0])):
        if a == values[0][i]:
            temp.append(values[1][i])
        else:
            list_moy.append(moy(temp)) 
            list_var.append(ecart_type(temp))
            abscisse.append(a)
            temp = []
            temp.append(values[1][i])
            a = values[0][i]
    list_moy.append(moy(temp)) 
    list_var.append(ecart_type(temp))
    abscisse.append(a)
    
    return(abscisse,list_moy,list_var)

            
def moyenne_glissante(pos,force):
    k = 7500
    pos_cut = pos[k:len(pos)-k-1]
    moy_force = []
    somme = 0
    for i in range(len(force)):
        if i < 2*k+1:
            somme += force[i]
        elif i == 2*k+1:
            somme += force[i]
            moy_force.append(somme/(2*k+1))
        else:
            somme += force[i]
            somme -= force[i-2*k-1]
            moy_force.append(somme/(2*k+1))
    return(pos_cut,moy_force)
    
    
def plot_curve(path_list,title,axe_list,commentaire):
    print(f"Titre :  {title}")
    print(f"Axes : {axe_list}")
    for path in path_list:
        print(path)
    print(commentaire)
    values = [[] for i in range(len(axe_list))]    
    for path in path_list:
        file = lk.File(path)
        path = path.split("\\")[-1]
        try:
            print(path+"        1")
            fd = file.fdcurves[list(file.fdcurves)[0]]
            force = fd.f.data
            distance = fd.d.data
            list_d,list_f = traitement(distance,force)    
        except:
            try:
                print(path+"        2")
                pos = file["Trap position"]["1X"].data #/!\ CHANGEMENT DE CODE
                F = file["Force HF"]["Force 2x"].data
                #pos = file["distance"]['distance 1'].data
                #F = file['Force LF']["Force 2x"]
                plt.plot(pos,F)
                plt.title("raw_data")
                plt.show()
                force = [-i for i in F]
                pos_cut,moy_force = moyenne_glissante(pos,force)
                plt.plot(pos_cut,moy_force)
                plt.title("moyenne glissante")
                plt.show()
                list_d,list_f = traitement(pos_cut,moy_force)
            except:
                print(path+"        3 aie")
                pos = file["Distance"]["Distance 1"].data
                F = file["Force LF"]["Force 2x"].data
                force = [-i for i in F]
                list_d,list_f = traitement(pos,force)
        for i in range(len(list_d)):
            try:   #.\Data\Training force xy\20220401-163938 Marker PP23 c5-v0.5-t1-d10-p10-x1.h5
                if path.split(" ")[4].split('.')[-1] == "h5":
                    param_list = (".").join(path.split(" ")[4].split('.')[:len(path.split(" ")[4].split('.'))-1]).split("-")     
                else:
                    param_list = path.split(" ")[4].split("-")
            except:
                if path.split(" ")[3].split('.')[-1] == "h5":
                    param_list = (".").join(path.split(" ")[3].split('.')[:len(path.split(" ")[3].split('.'))-1]).split("-")     
                else:
                    param_list = path.split(" ")[3].split("-")
            print(f"param_list = {param_list}")
            for param in param_list:
                if param[0] == 'v':
                    v = float(param[1:])
                if param[0] == 't':
                    t = float(param[1:])  
                if param[0] == 'd':
                    d = float(param[1:])
                if param[0] == "x":
                    x = float(param[1:])
            plt.plot(list_d[i],list_f[i])
            plt.show()
            aire_hyst = Polygon([[list_d[i][j], list_f[i][j]] for j in range(len(list_d[i]))]).area   
            #Lp,Lc,fit=fit_WLC(list_f[i],list_d[i],f"v : {v}, t : {t}, c : {(i+1)}")
            try:
                pass         #,St,d_offset,f_offset               print(True)
            except:
                print("fit impossible")
            for k in range(len(axe_list)):
                axe = axe_list[k]
                if axe == "v":
                    values[k].append(v)
                elif axe == "t":
                    values[k].append(t)
                elif axe == "d":
                    values[k].append(d)
                elif axe == "force":
                    values[k].append(list_f[i])
                elif axe == "distance":
                    values[k].append(list_d[i])
                elif axe == "Lc":
                    values[k].append(Lc)
                elif axe == "Lp":
                    values[k].append(Lp)
                elif axe == "aire_hyst":
                    values[k].append(aire_hyst)
                elif axe == 'x':
                    values[k].append(math.log(x,10))
                elif axe == 'fit':
                    values[k].append(fit) #(Lp,Lc,St,d_offset,f_offset))
    print("Commentaire : " + commentaire)        
    if axe_list[0]!="distance":
        if len(axe_list) == 3:
            fig, ax1 = plt.subplots(figsize=(16, 12))
            ax1.plot(values[0],values[1],".",color = "r", label = axe_list[1],markersize = 20)
            ax2 = ax1.twinx()
            ax2.plot(values[0],values[2],".",color = "b", label = axe_list[2],markersize = 20)
            ax1.set_ylabel(axe_list[1],fontsize = 20, fontweight='bold')
            ax2.set_ylabel(axe_list[2],fontsize = 20, fontweight='bold')
            ax1.set_xlabel(axe_list[0],fontsize = 20, fontweight='bold')
            plt.title(title, fontsize = 20, fontweight='bold')
            fig.tight_layout()
            fig.legend(fontsize = 20)
            plt.show()
        elif axe_list[0] == "":
            moy1 = moy(values[1][0:4])
            moy2 = moy(values[1][5:9])
            var1 = ecart_type(values[1][0:4])
            var2 = ecart_type(values[1][5:9])
            fig = plt.figure(1, figsize=(16, 12))
            plt.errorbar([0,1], [moy1,moy2], yerr = [var1,var2],fmt = 'none', capsize = 10, ecolor = 'red', elinewidth = 2, capthick = 8)
            plt.plot([0,1],[moy1,moy2],".",color = "b", label = "Moyenne",markersize = 20)
            plt.ylabel(axe_list[1],fontsize = 20, fontweight='bold')
            #plt.xlabel("Fatigue", fontsize = 20, fontweight='bold')
            plt.ylim([0,max(values[1])+1]) #min(values[1])-1
            plt.xlim([-0.5,1.5])
            plt.xticks([0,1],["Pas fatigue","Fatigue"],fontsize = 20, fontweight='bold')
            plt.yticks(fontsize = 20)
            plt.xticks(fontsize = 20)
            plt.title(title, fontsize = 20, fontweight='bold')
            fig.legend(fontsize = 20)
            plt.show()
        #elif axe_list[0] == "":
            #fig = plt.figure(1, figsize=(16, 12))
            #plt.plot([0 for i in range(5)]+[1 for i in range(4)],values[1],".",color = "r", label = axe_list[1],markersize = 20)
            #plt.ylabel(axe_list[1],fontsize = 20, fontweight='bold')
            #plt.xlabel("Fatigue", fontsize = 20, fontweight='bold')
            #plt.yticks(fontsize = 20)
            #plt.xticks(fontsize = 20)
            #plt.title(title, fontsize = 20, fontweight='bold')
            #fig.legend(fontsize = 20)
            #plt.show()
        #elif len(axe_list) == 2:
            #fig = plt.figure(1, figsize=(16, 12))
            #plt.plot(values[0],values[1],".",color = "b", label = axe_list[1],markersize = 20)
            #plt.ylabel(axe_list[1],fontsize = 20, fontweight='bold')
            #plt.xlabel(axe_list[0], fontsize = 20, fontweight='bold')
            #plt.yticks(fontsize = 20)
            #plt.xticks(fontsize = 20)
            #plt.title(title, fontsize = 20, fontweight='bold')
            #fig.legend(fontsize = 20)
        elif len(axe_list) == 2:
            print(values)
            abscisse,list_moy,list_var = Moy_and_Var(values)
            abscisse,list_moy,list_var = tri(abscisse,list_moy,list_var)
            #absi,a,b = fiit(abscisse,list_moy)
            fig= plt.figure(1, figsize=(16, 12))
            ax = fig.add_subplot(111)
            #plt.text(0.5,0.5,f" a = {a[0][0]}\n b = {b[0]}",fontsize = 20,horizontalalignment='center',verticalalignment='center', transform = ax.transAxes)
            #plt.plot(absi,(a*absi+b)[0],color='grey')
            #plt.plot([0,0.5],[56,56],'grey')
            plt.errorbar(abscisse, list_moy, yerr = list_var,fmt = 'none', capsize = 10, ecolor = 'red', elinewidth = 2, capthick = 8)            
            plt.plot(abscisse,list_moy,".",color = "b", label = "Moyenne",markersize = 20)
            plt.ylabel(axe_list[1],fontsize = 20, fontweight='bold')
            plt.xlabel(axe_list[0], fontsize = 20, fontweight='bold')
            plt.ylim([0,max(values[1])])
            plt.yticks(fontsize = 20)
            plt.xticks(fontsize = 20)
            plt.title(title, fontsize = 20, fontweight='bold')
            fig.legend(fontsize = 20)
            plt.show()
            if input("Save ? [y,n] : ")=='y':    #/!\---------------------------------------------------------------
                x = "5*10$^{-3}$mol/L"
                file = open("THE ENDING OF IT ALL.txt",'a')
                for compteur in range(len(abscisse)):
                    string = f'{x};{t};{abscisse[compteur]};{list_moy[compteur]};{list_var[compteur]}'
                    file.write(string+"\n")
                file.close()
    elif axe_list[0] == "distance" and axe_list[1] =='force' and len(axe_list) == 2:
        c=0
        fig = plt.figure(1, figsize=(16, 12))
        for i in range(len(path_list)):
            if path_list[i].split(" ")[4].split('.')[-1] == "h5":
                param = (".").join(path_list[i].split(" ")[4].split('.')[:len(path_list[i].split(" ")[4].split('.'))-1])
            else:
                param = path_list[i].split(" ")[4]
            for parametre in param.split('-'):
                if parametre[0] == 't':
                    t = int(parametre[1:]) 
            point = 100
            distance = values[0][c][:point]
            force = values[1][c][:point]
            qualibration = force[0] 
            force_qualibrée = [f - qualibration for f in force]
            plt.plot(distance,force_qualibrée,".-",linewidth = 3,color = color[i])
            plt.plot(distance,force_qualibrée,".",label = f"x : {x}",markersize = 10,color = color[i])
            
            for parametre in param.split('-'):
                if parametre[0] == 'c':
                    c += int(parametre[1:]) 
        plt.ylabel(axe_list[1],fontsize = 20, fontweight='bold')
        plt.xlabel(axe_list[0],fontsize = 20, fontweight='bold')
        plt.yticks(fontsize = 20)
        plt.xticks(fontsize = 20)
        plt.title(title, fontsize = 20, fontweight='bold')
        plt.legend(fontsize = 20)
        plt.show()    
    elif axe_list[0] == "distance" and axe_list[1] =='force' and axe_list[2] == 'fit':
        j = 0
        for i in range(len(path_list)):
            if path_list[i].split(" ")[4].split('.')[-1] == "h5":
                param = (".").join(path_list[i].split(" ")[4].split('.')[:len(path_list[i].split(" ")[4].split('.'))-1])
            else:
                param = path_list[i].split(" ")[4]
            for parametre in param.split('-'):
                if parametre[0] == 'c':
                    c = int(parametre[1:])
                if param[0] == 'v':
                    v = float(param[1:])
                if param[0] == 't':
                    t = float(param[1:]) 
            #plt.plot(values[0][i],values[1][i],label = param)
            for k in range(c):
                distance2,force2 = cut(values[0][j],values[1][j])
                fit = values[2][j] #wlc_fit = [4.11/(4*Lp)*((1-d/Lc)**(-2)-1+4*d/Lc) for d in distance2]
                fit.plot()#data = f"v : {v}, t : {t}, c : {i}")
                j += 1
        plt.legend()
        plt.show()
    

#f(d) = argmin[f](norm(DNA.Lc * (1 - (1/2)*sqrt(kT/(f*DNA.Lp)) + f/DNA.St)-(d - DNA.d_offset))) + DNA.f_offset
def moy(liste):
    return(sum(liste)/len(liste))
def ecart_type(liste):
    e = 0
    m = moy(liste)
    for i in range(len(liste)):
        e += (m - liste[i])**2
    e = e/len(liste)
    return(math.sqrt(e))
    
            
def save_a_plot(path_list,title,axe_list,commentaire):
    file = open("plot.txt",'a')
    string = ''
    for i in range(len(path_list)):
        if i != len(path_list)-1:
            string += path_list[i] + ';'
        else:
            string += path_list[i] + '    '
    string+=title + '    '
    for i in range(len(axe_list)):
        if i != len(axe_list)-1:
            string += axe_list[i] + ';'
        else:
            string += axe_list[i] + '    '
    string+=commentaire+'\n'
    file.write(string)
    file.close()
    
def open_and_show(choix = -1):
    if choix == -1:
        file = open("plot.txt",'r')
        line_list = file.readlines()
        file.close()
        for i in range(len(line_list)):
            print("("+str(i)+") "+line_list[i].split('    ')[1])
        curve = int(input("Curve number : "))
    else:
        file = open("plot.txt",'r')
        line_list = file.readlines()
        file.close()
        curve = choix
    path_list = line_list[curve].split('    ')[0].split(';')
    title = line_list[curve].split('    ')[1]
    axe_list = line_list[curve].split('    ')[2].split(';')
    commentaire = line_list[curve].split('    ')[3]
    plot_curve(path_list,title,axe_list,commentaire)
    
print("Okay")
    


# In[31]:


import os
list_h5 = []
for root, dirs, files in os.walk(".", topdown=True):
    print(dirs)
    for name in files:
        if name.split(".")[-1] == "h5":
            list_h5 += [os.path.join(root, name) ]
print(list_h5)
i = 0
file = list_h5[0].split("\\")[2]
print('-'+file+'------------------------------------')
for path in list_h5:
    if path.split("\\")[2] != file:
        print('-'+path.split("\\")[2]+'------------------------------------')
        file = path.split("\\")[2]
    print("("+str(i)+") "+(' ').join(path.split('\\')[-1].split(' ')[1:]))
    i+=1
Nb = int(input("Fichier numéro : "))
path = list_h5[Nb]

file = lk.File(path)



print(file)


# In[22]:


#MAKE A PLOT

def make_a_plot(list_h5):
    curve_list = [552,553,554,555,556,557,558,559,560,561,562]
    title = "Lp(x) : t = 0, v = 1"
    axe_list = ["v",'aire_hyst'] #[distance,force,v,t,Lc,Lp,aire_hyst,fit,""]
    commentaire = "Mesure de l'aire de l'hystèresis pour t = 0s et x = 1"
        
    path_list = []
    for curve in curve_list:
        path_list.append(list_h5[curve])

    plot_curve(path_list,title,axe_list,commentaire)
    
    if input("good ? [y,n] : ")=="y":
        save_a_plot(path_list,title,axe_list,commentaire)
        
make_a_plot(list_h5)


# In[ ]:





# In[ ]:





# In[ ]:





# In[159]:


def moyenne_5curve(pos,force):
    N = min([len(p) for p in pos])
    curve_moy = []
    force_moy = []
    for i in range(N):
        curve_moy.append(moy([p[i] for p in pos]))
        force_moy.append(moy([f[i] for f in force]))
    force_zero = []
    for i in range(N):
        force_zero.append(force_moy[i]-force_moy[0])
    return(curve_moy, force_zero)





def plot_5curves(Nb):

    path = list_h5[Nb]


    file = lk.File(path)

    pos = file["Distance"]["Distance 1"].data #/!\ CHANGEMENT DE CODE
    F = file["Force LF"]["Force 2x"].data
    force = [-f for f in F]

    pos_trait, force_trait = traitement(pos,force)
    pos_moy, force_moy = moyenne_5curve(pos_trait,force_trait)
    plt.plot(pos_moy,force_moy,'.', markersize = 20, label = label_list[Nb_list.index(Nb)])
Nb_list = [645,237,279]           
label_list = ["[Na$^+$] = 5 Mm", "[Na$^+$] = 0.5 Mm","[Na$^+$] < 0.05 Mm"]             
fig = plt.figure(1, figsize=(16, 12))
for Nb in [645,237,279]:
    plot_5curves(Nb)
plt.ylabel("Force (pN)",fontsize = 40)
plt.xlabel("Elongation ($\mu$m)", fontsize = 40)
plt.yticks(fontsize = 40)
plt.xticks(fontsize = 40)
plt.legend(fontsize = 40)
plt.savefig("Final Courbe Lp.png",dpi = "figure",format="png",transparent = True)
plt.show()
  


# In[145]:


#SHOW A PLOT

open_and_show(choix = -1)


# In[8]:


#TRAITEMENT MARKER AUTOMATIQUE

list_courbe_relax = [i for i in range(491,511)]#+[i for i in range(510,515)]#+[i for i in range(382,385)]


d_list = [-7.5,-12.5,-17.5,-20]
mLF_list = [7.5,12.5,17.5,20]
f = 78125

l = 0

for Nb in list_courbe_relax:
    path = list_h5[Nb]
    file = lk.File(path)

    def erfunc(x,mLF,a,b,d):
        return(mLF*scipy.special.erf((x-a)/(b*np.sqrt(2)))+d)


    force = file["Force HF"]["Force 2x"].data
    times = file["Force HF"]["Force 2x"].timestamps

    #seuil = moy(force[:642750]) # Marker 2
    #seuilseuil = moy(force[643250:]) # Marker 2

    N1 = 100
    N2 = 100
    seuil = moy(force[:N1])  # Marker 1
    seuilseuil = moy(force[N2:])  # Marker 1

    seuil1 = [seuil for i in range(len(force))] 
    seuil2 = [seuilseuil for i in range(len(force))]

    
    #fig = plt.figure(1, figsize=(16, 12))
    #plt.plot(force)
    #plt.show()
    print([mLF_list,[int(l//5)],len(force)/2,100,d_list[int(l//5)]])
    #params, extras = scipy.optimize.curve_fit(erfunc, np.array([i for i in range(len(force))]), force, p0=[25,643080,120,-25]) # Marker 2
    try:
        params, extras = scipy.optimize.curve_fit(erfunc, np.array([i for i in range(len(force))]), force, p0=[mLF_list[int(l//5)],len(force)/2,100,d_list[int(l//5)]]) # Marker 1
        print(params)
        [mLF,a,b,d] = params

    except:

        print("error1")
    print(f"l = {l}")
    l+=1
    fig = plt.figure(1, figsize=(16, 12))

    #plt.xlim([0.642*1e6,0.644*1e6])# Marker 2



    plt.ylabel("Force",fontsize = 20, fontweight='bold')
    plt.xlabel("Nb de point", fontsize = 20, fontweight='bold')
    plt.yticks(fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.title("Courbe de relaxation du brin d'ADN", fontsize = 20, fontweight='bold')

    plt.plot(force,".",linewidth = 5,label = "Force")

    try:
        plt.plot(erfunc(np.array([i for i in range(len(force))]),*params),label = "Erf")

        plt.plot([a,a],[d-mLF-5,10],'black')
        plt.plot([a-b*1.44*2,a-b*1.44*2],[d-mLF-5,10],'red')
        plt.plot([a+b*1.44*2,a+b*1.44*2],[d-mLF-5,10],'red')
        plt.plot([a-b*1.44,a-b*1.44],[d-mLF-5,10],'red')
        plt.plot([a+b*1.44,a+b*1.44],[d-mLF-5,10],'red')
        plt.plot([-100,len(force)],[d-mLF,d-mLF])
        plt.xlim([0,len(times)-1])
        fig.legend(fontsize = 20)
        plt.show()
        F = d-mLF
        tau  = b*1.44*4/f
        print(f"F = {F}")
        print(f"tau = {tau}")
    except:
        print("error2")
    print(f'Nb = {Nb}')
    save = 1
    if save:
        file = open("Marker x0.1 new.txt","a")
        file.write(f"{Nb};{a};{tau};{F};{path}\n")
        file.close()


# In[113]:


#TRAITEMENT MARKER AUTOMATIQUE NOUVELLE VERSION
from math import exp

#list_courbe_relax = [i for i in range(491,492)]#+[i for i in range(510,515)]#+[i for i in range(382,385)]
d_list = [-7.5,-12.5,-17.5,-20]
mLF_list = [7.5,12.5,17.5,20]
f = 78125
k_trap_ulrich = 2.3e-5
k_trap = 0.18e-3
eta = 10e-3
F0p = 40e-12
Lp = 56e-9
Lc = 16.5e-6
R = 4.98e-6
a = 2e-9
kBT= 1.38e-23*300
pi = 3.141592653
zeta_ADN = 2*pi*eta/math.log(Lc/a)
gamma_bead = 6*pi*eta*R
print(zeta_ADN*Lc,gamma_bead)

alpha = (gamma_bead+zeta_ADN*Lc)/gamma_bead
xi = (zeta_ADN*Lc)/gamma_bead/4*k_trap*Lc/(f0p**(3/2))*(kBT/Lp)**(1/2)
print(f"alpha = {alpha}, xi = {xi}")
print(math.log(alpha,10),math.log(xi,10))

beta = (xi/alpha)**0.56/((xi/alpha)**0.56+0.6*(alpha/xi)**0.56)*alpha/(1+alpha)+1/(1+alpha)
tau = (1+alpha)/xi+1.7*xi/alpha+0.16/(xi*alpha)
print(tau,beta)

t0 = (gamma_bead+zeta_ADN*Lc)/4*(kBT/Lp)**(1/2)*Lc/f0p**(3/2)
print(t0)

l = 0

def fp(t):
    return(beta*F0p/(1+9*(t/t0)**2)**(1/3)+(1-beta)*f0p*exp(-t/(t0*tau)))

t = [i*0.00001 for i in range(100000)]
force_theorique = [fp(time) for time in t]

plt.plot(t,force_theorique)
plt.show()
for Nb in list_courbe_relax:
    path = list_h5[Nb]
    file = lk.File(path)

    def erfunc(x,mLF,a,b,d):
        return(mLF*scipy.special.erf((x-a)/(b*np.sqrt(2)))+d)


    force = file["Force HF"]["Force 2x"].data
    times = file["Force HF"]["Force 2x"].timestamps
    F = [-f for f in force]
    #seuil = moy(force[:642750]) # Marker 2
    #seuilseuil = moy(force[643250:]) # Marker 2

 
    
    #fig = plt.figure(1, figsize=(16, 12))
    #plt.plot(force)
    #plt.show()
    print([mLF_list,[int(l//5)],len(force)/2,100,d_list[int(l//5)]])
    #params, extras = scipy.optimize.curve_fit(erfunc, np.array([i for i in range(len(force))]), force, p0=[25,643080,120,-25]) # Marker 2
    try:
        params, extras = scipy.optimize.curve_fit(erfunc, np.array([i for i in range(len(force))]), force, p0=[mLF_list[int(l//5)],len(force)/2,100,d_list[int(l//5)]]) # Marker 1
        print(params)
        [mLF,a,b,d] = params

    except:

        print("error1")
    print(f"l = {l}")
    l+=1
    fig = plt.figure(1, figsize=(16, 12))

    #plt.xlim([0.642*1e6,0.644*1e6])# Marker 2



    plt.ylabel("Force",fontsize = 20, fontweight='bold')
    plt.xlabel("Nb de point", fontsize = 20, fontweight='bold')
    plt.yticks(fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.title("Courbe de relaxation du brin d'ADN", fontsize = 20, fontweight='bold')

    plt.plot(force,".",linewidth = 5,label = "Force")

    try:
        plt.plot(erfunc(np.array([i for i in range(len(force))]),*params),label = "Erf")

        plt.plot([a,a],[d-mLF-5,10],'black')
        plt.plot([a-b*1.44*2,a-b*1.44*2],[d-mLF-5,10],'red')
        plt.plot([a+b*1.44*2,a+b*1.44*2],[d-mLF-5,10],'red')
        plt.plot([a-b*1.44,a-b*1.44],[d-mLF-5,10],'red')
        plt.plot([a+b*1.44,a+b*1.44],[d-mLF-5,10],'red')
        plt.plot([-100,len(force)],[d-mLF,d-mLF])
        plt.xlim([0,len(times)-1])
        fig.legend(fontsize = 20)
        plt.show()
        F = d-mLF
        tau  = b*1.44*4/f
        print(f"F = {F}")
        print(f"tau = {tau}")
    except:
        print("error2")
    print(f'Nb = {Nb}')
    save = 0
    if save:
        file = open("Marker x0.1 new.txt","a")
        file.write(f"{Nb};{a};{tau};{F};{path}\n")
        file.close()


# In[117]:


#PLOT TEMPS DE RELAXATION



import math

e_t = 1

l_o_p = 0
matrice_relax = []
matrice_force = []

def log_ou_pas(x):
    if l_o_p:
        return(math.log(x,10))
    else:
        return(x)
    

fig = plt.figure(1, figsize=(16, 12))

file = open("Marker x1 new.txt","r")
L = file.readlines()
file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))

list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))
    
    
matrice_force += list_moy_F
matrice_relax += list_moy_T
color_curve='red'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]

reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)
    
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 5 mM",color=color_curve)


#---------------------------------------------------------------------------------------------------------------------


file = open("Marker x0.1 new.txt","r")
L = file.readlines()
file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))

list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))
    

matrice_force += list_moy_F
matrice_relax += list_moy_T

color_curve='blue'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.5 mM",color=color_curve)



#---------------------------------------------------------------------------------------------------------------------
file = open("Marker x1 dans le flux.txt","r")
L = file.readlines()

file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))
   
    
list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):  #Nombre de jeu de données /!\
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))


matrice_force += list_moy_F
matrice_relax += list_moy_T
   
color_curve='orange'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.25 mM",color=color_curve)

#---------------------------------------------------------------------------------------------------------------------


file = open("Marker x0.05.txt","r")
L = file.readlines()

file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))
   
    
list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):  #Nombre de jeu de données /!\
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))

matrice_force += list_moy_F
matrice_relax += list_moy_T
    
color_curve='green'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)   
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.1 mM",color=color_curve)



#---------------------------------------------------------------------------------------------------------------------


file = open("Marker x0.01 bis sans lissage.txt","r")
L = file.readlines()

file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))
   
    
list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(8):  #Nombre de jeu de données /!\
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))

matrice_force += list_moy_F
matrice_relax += list_moy_T
  
color_curve='purple'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)  
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.05 mM",color=color_curve)

#---------------------------------------------------------------------------------------------------------------------



[absi,a,b] = fiit(matrice_force,matrice_relax)
a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+11) for i in range(44)]
reg_time = [a*r_f+b for r_f in reg_force]
plt.plot(reg_force,reg_time,'--',color='grey',linewidth = 5,label = "Fit : a = 7.58e-2 ms/pN")




#f0p = 40e-12
#list_f0p = [15e-12,20e-12,25e-12,30e-12,35e-12,40e-12,45e-12,50e-12]
#list_f0p_l_o_p = [log_ou_pas(15),log_ou_pas(20),log_ou_pas(25),log_ou_pas(30),log_ou_pas(35),log_ou_pas(40),log_ou_pas(45),log_ou_pas(50)]
#print(len(list_f0p),len(list_f0p_l_o_p))
#fig = plt.figure(1,figsize = (16, 12))
#for coeff in [0.5,0.4,0.3,0.2,0.1]:
#    list_t_relax = []
#    for f0p in list_f0p:
#        list_t_relax.append(resoudre_Fx_v(f0p,coeff))
#    
#    [absi,a,b] = fiit(list_f0p_l_o_p,list_t_relax)
#    print(f'a = {a[0][0]}, b = {b[0]}')##
#
#    a = a[0][0]
#    b = b[0]
#    reg_force = list_f0p_l_o_p
#    reg_time = [a*r_f+b for r_f in reg_force]
#    if l_o_p:
#        plt.plot(reg_force,reg_time,'-.',color='grey')  
#    
#    plt.plot(list_f0p_l_o_p,list_t_relax,'.',label = str(coeff))







plt.ylabel("temps de relaxation (ms)",fontsize = 20, fontweight='bold')
plt.xlabel("Force de maintient (pN)", fontsize = 20, fontweight='bold')
#plt.ylim([0,max(values[1])])
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 20)
plt.ylim([0.0075,0.0155])
plt.xlim([10,65])
plt.yticks([0.008,0.009,0.010,0.011,0.012,0.013,0.014,0.015],["8","9","10","11","12","13","14","15"])
plt.legend(loc= "lower left",fontsize = 18)
plt.savefig("final Temps de relax.png",dpi = "figure",format="png",transparent = True)
plt.show()
    


# In[24]:


#PLOT TEMPS DE RELAXATION



import math

e_t = 0

l_o_p = 1
matrice_relax = []
matrice_force = []

def log_ou_pas(x):
    if l_o_p:
        return(math.log(x,10))
    else:
        return(x)
    

fig = plt.figure(1, figsize=(16, 12))

file = open("Marker x1 new.txt","r")
L = file.readlines()
file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))

list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))
    
    
matrice_force += list_moy_F
matrice_relax += list_moy_T
color_curve='red'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]

reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)
    
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 5 mM",color=color_curve)


#---------------------------------------------------------------------------------------------------------------------


file = open("Marker x0.1 new.txt","r")
L = file.readlines()
file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))

list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))
    

matrice_force += list_moy_F
matrice_relax += list_moy_T

color_curve='blue'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.5 mM",color=color_curve)



#---------------------------------------------------------------------------------------------------------------------
file = open("Marker x1 dans le flux.txt","r")
L = file.readlines()

file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))
   
    
list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):  #Nombre de jeu de données /!\
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))


matrice_force += list_moy_F
matrice_relax += list_moy_T
   
color_curve='orange'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.25 mM",color=color_curve)

#---------------------------------------------------------------------------------------------------------------------


file = open("Marker x0.05.txt","r")
L = file.readlines()

file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))
   
    
list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(4):  #Nombre de jeu de données /!\
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))

matrice_force += list_moy_F
matrice_relax += list_moy_T
    
color_curve='green'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)   
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.1 mM",color=color_curve)



#---------------------------------------------------------------------------------------------------------------------


file = open("Marker x0.01 bis sans lissage.txt","r")
L = file.readlines()

file.close()
list_F = []
list_T = []
for i in range(len(L)):
    list_F.append(float(L[i].split(";")[3]))
    list_T.append(float(L[i].split(";")[2]))
   
    
list_moy_F = []
list_moy_T = []
list_var_F = []
list_var_T = []

for i in range(8):  #Nombre de jeu de données /!\
    list_moy_F.append(log_ou_pas(abs(moy(list_F[i*5:(i+1)*5]))))
    list_moy_T.append(log_ou_pas(moy(list_T[i*5:(i+1)*5])))
    list_var_F.append(ecart_type(list_F[i*5:(i+1)*5]))
    list_var_T.append(log_ou_pas(ecart_type(list_T[i*5:(i+1)*5])))

matrice_force += list_moy_F
matrice_relax += list_moy_T
  
color_curve='purple'    
[absi,a,b] = fiit(list_moy_F,list_moy_T)
print(f'a = {a}, b = {b}')

a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]

if e_t : 
    plt.errorbar([abs(i) for i in list_moy_F], list_moy_T, yerr = list_var_T,fmt = 'none', capsize = 10, ecolor = color_curve, elinewidth = 2, capthick = 8)
else:
    plt.plot(reg_force,reg_time,'-.',color=color_curve)  
plt.plot([abs(i) for i in list_moy_F],list_moy_T,".",markersize=20,label="[Na$^+$] = 0.05 mM",color=color_curve)

#---------------------------------------------------------------------------------------------------------------------



[absi,a,b] = fiit(matrice_force,matrice_relax)
a = a[0][0]
b = b[0]
reg_force = [log_ou_pas(i+10) for i in range(45)]
reg_time = [a*r_f+b for r_f in reg_force]
#plt.plot(reg_force,reg_time,'--',color='grey',linewidth = 5,label = "Fit : a = 7.58e-2 ms/pN")




#f0p = 40e-12
#list_f0p = [15e-12,20e-12,25e-12,30e-12,35e-12,40e-12,45e-12,50e-12]
#list_f0p_l_o_p = [log_ou_pas(15),log_ou_pas(20),log_ou_pas(25),log_ou_pas(30),log_ou_pas(35),log_ou_pas(40),log_ou_pas(45),log_ou_pas(50)]
#print(len(list_f0p),len(list_f0p_l_o_p))
#fig = plt.figure(1,figsize = (16, 12))
#for coeff in [0.5,0.4,0.3,0.2,0.1]:
#    list_t_relax = []
#    for f0p in list_f0p:
#        list_t_relax.append(resoudre_Fx_v(f0p,coeff))
#    
#    [absi,a,b] = fiit(list_f0p_l_o_p,list_t_relax)
#    print(f'a = {a[0][0]}, b = {b[0]}')##
#
#    a = a[0][0]
#    b = b[0]
#    reg_force = list_f0p_l_o_p
#    reg_time = [a*r_f+b for r_f in reg_force]
#    if l_o_p:
#        plt.plot(reg_force,reg_time,'-.',color='grey')  
#    
#    plt.plot(list_f0p_l_o_p,list_t_relax,'.',label = str(coeff))







plt.ylabel("temps de relaxation (ms)",fontsize = 20, fontweight='bold')
plt.xlabel("Force de maintient (pN)", fontsize = 20, fontweight='bold')
#plt.ylim([0,max(values[1])])
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 20)
#plt.ylim([0.008,0.014])
#plt.yticks([0.008,0.009,0.010,0.011,0.012,0.013,0.014],["8","9","10","11","12","13","14"])
plt.legend(loc= "lower left",fontsize = 20)
plt.savefig("final Temps de relax.png",dpi = "figure",format="png",transparent = True)
plt.show()
    


# In[109]:


#FIT UNE RELAXATION PARTICULIERE

import scipy.special
import scipy.optimize
import numpy as np
def erfunc(x,mLF,a,b,d):
    return(mLF*scipy.special.erf((x-a)/(b*np.sqrt(2)))+d)

f = 78125
force = file["Force HF"]["Force 2x"].data
times = file["Force HF"]["Force 2x"].timestamps
temps = [i/f for i in range(len(force))]

plt.plot([t-times[0] for t in times],temps)
plt.show()


#params, extras = scipy.optimize.curve_fit(erfunc, np.array([i for i in range(len(force))]), force, p0=[25,643080,120,-25]) # Marker 2
try:
    params, extras = scipy.optimize.curve_fit(erfunc, np.array(temps), force, p0=[25,250,120,-25]) # Marker 1
    print(params)
    [mLF,a,b,d] = params

except:
    print("error1")

fig = plt.figure(1, figsize=(16, 12))
ax = plt.gca()
#plt.xlim([0.642*1e6,0.644*1e6])# Marker 2



plt.ylabel("Force (pN)",fontsize = 40)
plt.xlabel("Temps (ms)", fontsize = 40)
plt.yticks(fontsize = 40)
plt.xticks(fontsize = 40)

ax.plot(temps,force,".",linewidth = 5,label = "Données expérimentales",markersize = 20)
ax.plot([-100,len(force)],[d-mLF,d-mLF],label = r'$F_{imposee}$',color = "green",linewidth = 5)
ax.plot(temps,erfunc(np.array(temps),*params),label = r"Fit : $\frac{F_{imposee}}{2}[erf(\frac{x-a}{b*\sqrt{2}}$)-1]",color= 'orange',linewidth = 5)

ax.plot([a,a],[d-mLF-5,10],'black',label = 'Fit : a = 18 ms',linewidth = 5)
ax.plot([a-b*1.44*2,a-b*1.44*2],[d-mLF-5,10],'red',label = "Fit : t$_{relax}$ = 4*b*$\sqrt{2}$ = 9.8 ms",linewidth = 5)
ax.plot([a+b*1.44*2,a+b*1.44*2],[d-mLF-5,10],'red',linewidth = 5)
plt.xticks([0,0.005,0.010,0.015,0.020,0.025,0.030],[0,5,10,15,20,25,30])
plt.xlim([0,0.03])
plt.legend(loc="upper left",fontsize = 30)

plt.savefig("Final temps de relax une courbe.png",dpi = "figure",format="png",transparent = True)

plt.show()
F = d-mLF
tau  = b*1.44*4
print(f"F = {F}")
print(f"tau = {tau}")


# In[71]:


#SAVE UNE RELAXATION PARTICULIERE

file = open("Marker x0.01.txt","a")
file.write(f"{Nb};{a};{tau*3/4};{F};{path}\n")
file.close()


# In[41]:


#PLOT COURBE HYST DEPUIS UN FICHIER

from colour import Color
red = Color("orange")
colors = list(red.range_to(Color("blue"),4))

file = open("data_hysteresis - Copie (3).txt",'r')
lines = file.readlines()
file.close()

Matrice_donnée = []
for line in lines:
    data = line.split(";")
    Matrice_donnée.append([data[0],float(data[1]),float(data[2]),float(data[3]),float(data[4])])

list_concentration = []
for data in Matrice_donnée:
    if not(data[0] in list_concentration):
        list_concentration.append(data[0])
print(list_concentration)
list_concentration_label=["[Na$^+$] = 5 mM","[Na$^+$] = 2.5 mM","[Na$^+$] = 1 mM","[Na$^+$] = 0.5 mM"]
list_a_label = ["Fit : a = 0.141","Fit : a = 0.156","Fit : a = 0.157","Fit : a = 0.173"]
t = 0

list_abs = []
list_ord = []
list_error = []

fig= plt.figure(1, figsize=(16, 12))

l_o_p = 1

def log_ou_pas(x):
    if l_o_p:
        return(math.log(x,10))
    else:
        return(x)

color_list = ["blue","orange","green","purple"]
i = 0
for x in list_concentration:
    abscisse = []
    ordonnée = []
    error_bar = []
    for data in Matrice_donnée:
        if data[0] == x:
            if data[1] == t:
                abscisse.append(log_ou_pas(data[2]))
                ordonnée.append(log_ou_pas(data[3]))
                error_bar.append(data[4])
    list_abs.append(abscisse)
    list_ord.append(ordonnée)
    list_error.append(error_bar)
    (absi,a,b) = fiit(abscisse,ordonnée)
    print(f"a = {a[0][0]},b = {b[0]}")
    if not(l_o_p):
        plt.errorbar(abscisse, ordonnée, yerr = error_bar,fmt = 'none', capsize = 10, ecolor = color_list[list_concentration.index(x)], elinewidth = 2, capthick = 8)            
    else:
        plt.plot(absi,(a*absi+b)[0],"--",color=color_list[list_concentration.index(x)],label = list_a_label[i],linewidth =5)
    plt.plot(abscisse,ordonnée,".",color = color_list[list_concentration.index(x)],markersize = 40)
    i+=1
plt.ylabel("$Aire~de~l'hystérésis~(log)$",fontsize = 30, fontweight='bold')
plt.xlabel("$Vitesse~(log)$", fontsize = 30, fontweight='bold')

plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.legend(fontsize = 30)
plt.savefig("final aire_hyst(x) log.png",dpi = "figure",format="png",transparent = True)
plt.show()


            
    



# In[48]:


#PLOT COURBE HYST DEPUIS UN FICHIER

from colour import Color
red = Color("orange")
colors = list(red.range_to(Color("blue"),4))

file = open("data_hysteresis - Copie (3).txt",'r')
lines = file.readlines()
file.close()

Matrice_donnée = []
for line in lines:
    data = line.split(";")
    Matrice_donnée.append([data[0],float(data[1]),float(data[2]),float(data[3]),float(data[4])])

list_concentration = []
for data in Matrice_donnée:
    if not(data[0] in list_concentration):
        list_concentration.append(data[0])
print(list_concentration)
list_concentration_label=["[Na$^+$] = 5 mM","[Na$^+$] = 2.5 mM","[Na$^+$] = 1 mM","[Na$^+$] = 0.5 mM"]
t = 0

list_abs = []
list_ord = []
list_error = []

fig= plt.figure(1, figsize=(16, 12))

l_o_p = 0

def log_ou_pas(x):
    if l_o_p:
        return(math.log(x,10))
    else:
        return(x)
j = 0
color_list = ["blue","orange","green","purple"]
for x in list_concentration:
    abscisse = []
    ordonnée = []
    error_bar = []
    for data in Matrice_donnée:
        if data[0] == x:
            if data[1] == t:
                abscisse.append(log_ou_pas(data[2]))
                ordonnée.append(log_ou_pas(data[3]))
                error_bar.append(data[4])
    list_abs.append(abscisse)
    list_ord.append(ordonnée)
    list_error.append(error_bar)
    (absi,a,b) = fiit(abscisse,ordonnée)
    print(f"a = {a[0][0]},b = {b[0]}")
    if not(l_o_p):
        plt.errorbar(abscisse, ordonnée, yerr = error_bar,fmt = 'none', capsize = 3, ecolor = 'grey', elinewidth = 2, capthick = 8)            
    else:
        plt.plot(absi,(a*absi+b)[0],color=color_list[list_concentration.index(x)],linewidth = 5)
    plt.plot(abscisse,ordonnée,".",color = color_list[list_concentration.index(x)], label = list_concentration_label[j],markersize = 20)
    j+=1
    #print(j,list_concentration_label[j] )
    
plt.ylabel("$Aire~de~l'hystérésis~(pN.\mu m)$",fontsize = 30, fontweight='bold')
plt.xlabel("$Vitesse~(\mu m/s)$", fontsize = 30, fontweight='bold')

plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.legend(loc = "lower right",fontsize = 25)
plt.xlim([-5,55]) 
plt.ylim([5,80]) 
plt.savefig("final aire_hyst(x).png",dpi = "figure",format="png",transparent = True)

plt.show()


            
    



# In[8]:


#PRINT FICHIER EN LOG

x = Matrice_donnée[0][0]
l_o_p = 1
list_t = []
for data in Matrice_donnée:
    if not(data[1] in list_t):
        list_t.append(data[1])
print(list_t)

list_abs = []
list_ord = []
list_error = []

fig= plt.figure(1, figsize=(16, 12))

color_list = ["blue","orange","green","purple"]
for t in list_t:
    abscisse = []
    ordonnée = []
    error_bar = []
    for data in Matrice_donnée:
        if data[1] == t:
            if data[0] == x:
                abscisse.append(log_ou_pas(data[2]))
                ordonnée.append(log_ou_pas(data[3]))
                error_bar.append(data[4])
    list_abs.append(abscisse)
    list_ord.append(ordonnée)
    list_error.append(error_bar)
    (absi,a,b) = fiit(abscisse,ordonnée)
    print(f"a = {a[0][0]},b = {b[0]}")
    if not(l_o_p):
        plt.errorbar(abscisse, ordonnée, yerr = error_bar,fmt = 'none', capsize = 10, ecolor = colors[list_t.index(t)].hex, elinewidth = 2, capthick = 8)            
    else:
        plt.plot(absi,(a*absi+b)[0],'--',color= colors[list_t.index(t)].hex,linewidth = 5, label = f"Fit : a = "+"{:0.3f}".format(a[0][0]))    
    plt.plot(abscisse,ordonnée,".",color = colors[list_t.index(t)].hex,markersize = 40)

plt.ylabel("$Aire~de~l'hystérésis~(log)$",fontsize = 40, fontweight='bold')
plt.xlabel("$Vitesse~(log)$", fontsize = 40, fontweight='bold')

plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.legend(fontsize = 40)
plt.savefig("final aire_hyst(t) log.png",dpi = "figure",format="png",transparent = True)
plt.show()


# In[135]:


x = Matrice_donnée[0][0]
l_o_p = 0
list_t = []
for data in Matrice_donnée:
    if not(data[1] in list_t):
        list_t.append(data[1])
print(list_t)

list_abs = []
list_ord = []
list_error = []

fig= plt.figure(1, figsize=(16, 12))

for t in list_t:
    abscisse = []
    ordonnée = []
    error_bar = []
    for data in Matrice_donnée:
        if data[1] == t:
            if data[0] == x:
                abscisse.append(log_ou_pas(data[2]))
                ordonnée.append(log_ou_pas(data[3]))
                error_bar.append(data[4])
    list_abs.append(abscisse)
    list_ord.append(ordonnée)
    list_error.append(error_bar)

for i in range(len(list_abs)):
    
    plt.errorbar(list_abs[i], list_ord[i], yerr = list_error[i],fmt = 'none', capsize = 3, ecolor = "grey", elinewidth = 2, capthick = 8)            
    
for i in range(len(list_abs)):
    plt.plot(list_abs[i],list_ord[i],".",color = colors[i].hex, label = f"temps de relaxation = {list_t[i]}s",markersize = 20)

plt.ylabel("$Aire~de~l'hystérésis~(pN.\mu m)$",fontsize = 30, fontweight='bold')
plt.xlabel("$Vitesse~(\mu m/s)$", fontsize = 30, fontweight='bold')

plt.yticks(fontsize = 30)
plt.xticks(fontsize = 30)
plt.legend(loc = "lower right",fontsize = 25)
plt.ylim([5,80])
plt.savefig("final aire_hyst(t).png",dpi = "figure",format="png",transparent = True) 
plt.show()


# In[123]:



def WLC_function(x,Lc,Lp,f_offset):
    #return(Lc*(1-1/2*(4.11/(F+offset)/Lp)**(1/2)+(F+offset)/Sm))
    return(4.11/Lp*(1/4*(1-x/Lc)**(-2)-1/4+x/Lc)+f_offset)

def fit_WLC(force,distance,path):
    d,f = cut(distance,force)
    print(len(d))
    min_f = min(f)
    j = f.index(min_f)

    force_cut = f
    distance_cut = d


    #for i in range(len(f)):
    #    if i<j:
     #       force_cut.append(f[j]-(f[i]-f[j]))
      #      distance_cut.append(d[i])
       # else:
        #    force_cut.append(f[i])
         #   distance_cut.append(d[i])
    
    print(len(distance_cut))
    f = []
    d = []
    for i in range(len(force_cut)):
        if force_cut[i]<10 and distance_cut[i]>12:
            f.append(force_cut[i])
            d.append(distance_cut[i])
    
    print(len(d))
    params, extras = scipy.optimize.curve_fit(WLC_function, np.array(d),np.array(f),p0=[16.5,56,-2]) # Marker 1
    fig  = plt.figure(1,figsize = (16, 12))
    plt.plot(d,f,'.',markersize = 20,label = "Données expérimentales")
    plt.plot(d,WLC_function(np.array(d),*params),'-.',label = "Modèle WLC Lp = 36",linewidth = 5)
    plt.plot(d,WLC_function(np.array(d),16.5,56,-0.43521548214540084),'--',label = "Modèle WLC Lp = 56",linewidth = 5)
    plt.ylabel("$Force~(pN)$",fontsize = 40, fontweight='bold')
    plt.xlabel("$Distance~(\mu m)$", fontsize = 40, fontweight='bold')
    plt.yticks(fontsize = 30)
    plt.xticks(fontsize = 30)
    plt.legend(fontsize = 40)
    plt.savefig("final Lp echec de fit.png",dpi = "figure",format="png",transparent = True)
    plt.show()
    print(f"Lc = {float(params[0])} f_offset = {float(params[2])}, Lp = {float(params[1])}")
    diff_fit = []
    for i in range(len(d)):
        diff_fit.append(-(WLC_function(np.array(d[i]),*params)-WLC_function(np.array(d[i]),16.5,56,-0.43521548214540084)))
    plt.plot(d,diff_fit)
    plt.show()
    return(float(params[1]),float(params[0]),float(params[2])) #
    
  

file = open("Data_lp.txt",'r')
data = file.readlines()
file.close()

force = []
distance = []

for line in data:
    distance.append(float(line.split(';')[0]))
    force.append(float(line.split(';')[1]))

fit_WLC(force,distance,path)



# In[2]:


#FIT WLC PAR LUMICKS

#d,f = cut(distance,force)
#model = lk.inverted_odijk("DNA").subtract_independent_offset() + lk.force_offset("DNA")
#fit = lk.FdFit(model)
#fit.add_data(path, f, d)
#fit['DNA/d_offset'].value = distance[0]
#fit['DNA/d_offset'].lower_bound = 0
#fit['DNA/d_offset'].upper_bound = distance[0]*1.5
#fit['DNA/f_offset'].value = force[0]
#fit['DNA/f_offset'].lower_bound = 0
#fit['DNA/f_offset'].upper_bound = force[0]*1.5
#fit['DNA/Lp'].value = 50
#fit['DNA/Lp'].lower_bound = 0
#fit['DNA/Lp'].upper_bound = 200
#fit['DNA/Lc'].value = 16
#fit['DNA/Lc'].lower_bound = 0
#fit['DNA/Lc'].upper_bound = 30
#fit['DNA/St'].value = 1500
#fit['DNA/St'].lower_bound = 0
#fit['DNA/St'].upper_bound = 2000
#fit.fit()
#print(model.get_formatted_equation_string(1))
#fit.plot()
#plt.ylabel("Force [pN]")#
#plt.xlabel("Distance [$\\mu$m]")
#plt.show()


# In[ ]:


#FIT WLC SCIPY

#import scipy.special
#import scipy.optimize
#import numpy as np
#def WLC_function(F,Lc,Lp,Sm,offset):
#    return(Lc*(1-1/2*(4.11/(F+offset)/Lp)**(1/2)+(F+offset)/Sm))


#fd = file.fdcurves[list(file.fdcurves)[0]]
#force = fd.f.data
#distance = fd.d.data
#        
#list_d,list_f = traitement(distance,force)
#
#d,f = cut(list_d[0],list_f[0])

#force_cut = []
#distance_cut = []

#min_f = min(f)
#j = f.index(min_f)

#for i in range(len(f)):
#    if i<j:
#        force_cut.append(f[j]-(f[i]-f[j]))
#        distance_cut.append(d[i])
#    else:
#        force_cut.append(f[i])
#        distance_cut.append(d[i])

#min_f = min(force_cut)

#temp = [force-min_f*1.2 for force in force_cut]

#force_cut = temp

#force_cut2 = []
#distance_cut2 = []

#for i in range(len(force_cut)):
#    if force_cut[i]>1:
#        force_cut2.append(force_cut[i])
#        distance_cut2.append(distance_cut[i])
        
        
        
#params, extras = scipy.optimize.curve_fit(WLC_function, np.array(force_cut),np.array(distance_cut),p0=[16.5,56,1500,0]) # Marker 1
#print(f"Lc = {float(params[0])} Lp = {params[1]} Sm = {params[2]} Offset = {params[3]}")
#fig  = plt.figure(1,figsize = (16, 12))
#plt.plot(distance_cut,force_cut)
#plt.plot(WLC_function(np.array(force_cut),*params),force_cut)
#plt.plot(WLC_function(np.array(f),*[16.46,56,0.46]),f,'red')

#plt.show()


# In[ ]:


#CALCUL AIRE HYSTERESIS

#from shapely.geometry import Polygon
#for j in range(len(list_d)):
#    print('('+str(j)+") "+str(Polygon([[list_d[j][i], list_f[j][i]] for i in range(len(list_d[j]))]).area))


# In[19]:


#PRINT UN PLOT PRECISEMENT 

#pathFJC = list_h5[88]

#fileFJC = lk.File(pathFJC)


#fd = file.fdcurves[list(fileFJC.fdcurves)[0]]
#forceFJC = fd.f.data
#distanceFJC = fd.d.data


pathWLC = list_h5[153]

fileWLC = lk.File(pathWLC)


distanceWLC = fileWLC["Trap position"]["1X"].data #/!\ CHANGEMENT DE CODE
forceWLC = fileWLC["Force HF"]["Force 2x"].data

plt.plot(distanceFJC,forceFJC)
plt.plot(distanceWLC,forceWLC)
plt.show()

Lp = 56
Lc = 16.5
Sm = 1500

def WLC_function_ext(F):
    return(Lc*(1-1/2*(4.11/F/Lp)**(1/2)+F/Sm))

def moyenne_glissante(pos,force):
    k = 4
    pos_cut = pos[k:len(pos)-k-1]
    moy_force = []
    somme = 0
    for i in range(len(force)):
        if i < 2*k+1:
            somme += force[i]
        elif i == 2*k+1:
            somme += force[i]
            moy_force.append(somme/(2*k+1))
        else:
            somme += force[i]
            somme -= force[i-2*k-1]
            moy_force.append(somme/(2*k+1))
    return(pos_cut,moy_force)
    

fig = plt.figure(1, figsize=(16, 12))
d,f = traitement(distance,force)

wlc= [WLC_function(dist,16.5,56,0) for dist in d[0][:200]]

a=200
b=232
c=344
e=447
d2=d[0]
f2=f[0]

min_f = min(f2)
j = f2.index(min_f)

force_ext =  [10+i*0.01 for i in range(7000)]
distance_ext = [WLC_function_ext(force_ext[i]) for i in range(len(force_ext))]


wlc= [WLC_function(dist,16.5,56,0) for dist in distance_cut]


# In[12]:



fig = plt.figure(1, figsize=(16, 12))
t=20
plt.plot([16.5,16.5],[-10,80],'--',linewidth = 8,color="darkgrey",label="Longueur de contour")
plt.plot([0],[0],'.',color = "#104E8B",markersize=t,label="Données inextensibles")
plt.plot([0],[0],'.',color = "#1E90FF",markersize=t,label="Données extensibles")
plt.plot([0],[0],'.',color = "#8B1A1A",markersize=t,label="Plateau de dénaturation")
plt.plot([0],[0],'.',color = "#EE2C2C",markersize=t,label="Hystérésis")
plt.plot([0],[0],color = "#E3CF57",linewidth = 5,label="WLC inextensible")
plt.plot([0],[0],color = "#D45500",linewidth = 5,label="WLC extensible")

plt.ylabel("$Force~(pN)$",fontsize = 30, fontweight='bold')
plt.xlabel("$Distance~(\mu m)$", fontsize = 30, fontweight='bold')
plt.yticks([i*5 for i in range(13)],[i*5 for i in range(13)],fontsize = 30)
plt.xticks([i for i in range(10,20)],[i for i in range(10,20)],fontsize = 30)
plt.legend(fontsize = 30)
plt.ylim([-5,65])
plt.xlim([10.5,19.5])
plt.savefig("final legend Essai en traction.png",dpi = "figure",format="png",transparent = True)
plt.show()


# In[28]:


#FIT WLC EXTENTIBLE AVEC LA FORCE ET SCIPY

#import scipy.special
#import scipy.optimize
#import numpy as np
#def WLC_function(F,Lp,Lc,Sm):
#    return(Lc*(1-1/2*(4.11/F/Lp)**(1/2)+F/Sm))


#fd = file.fdcurves[list(file.fdcurves)[0]]
#force = fd.f.data
#distance = fd.d.data
#list_d,list_f = traitement(distance,force)

#d,f = cut(list_d[0],list_f[0])

#params, extras = scipy.optimize.curve_fit(WLC_function, np.array(f),np.array(d),p0=[56e-3,16.5,1500]) # Marker 1
#print(params)
#fig  = plt.figure(1,figsize = (16, 12))
#plt.plot(d,f)
#plt.plot(WLC_function(np.array(f),*params),f)
#plt.plot(WLC_function(np.array(f),*[56e-3,16.5,1500]),f,'red')

#plt.show()


# In[29]:


#ANCIEN H5_LIST

#h5_list = []
#i = 0
#for path in os.listdir():
#    if path.split('.')[-1] == 'h5':
#        h5_list.append(path)
#        print("("+str(i)+") "+(' ').join(path.split(' ')[1:]))
#        i+=1
#Nb = int(input("Fichier numéro : "))
#path = h5_list[Nb]

#file = lk.File(path)
#print(file)


# In[162]:


import scipy.integrate 



def resoudre_Fx_v(f0p,coeff):
    eta = 1e-3
    R = 4.98e-6
    pi = 3.141592653
    Sm = 1500e-12
    Lp = 56e-9
    Lc = 16.5e-6
    f_relax = f0p*coeff
    masse_bille = 4/3*pi*R**3*997
    print(masse_bille)
    def WLC_function(F):
        return(Lc*(1-1/2*(1.38e-23*300/F/Lp)**(1/2)+F/Sm))
    def dFsurdt(F,t):
        return(-F/(6*pi*eta*R*(Lc/4*(1.38e-23*300/Lp)**(1/2)*F**(-3/2)+Lc/Sm)))
    t = [i*0.00005 for i in range(1000)]
    f_vect = scipy.integrate.odeint(dFsurdt,f0p,t)
    j=0
    while f_vect[j]>f_relax:
        j+=1
    t_relax = log_ou_pas(t[j])
    x_vect = []
    for i in range(len(f_vect)):
        x_vect.append(WLC_function(f_vect[i]))
    #fig = plt.figure(1,figsize = (16, 12))
    #plt.plot(t,x_vect)
    #plt.show()
    v_vect = []
    t_vitesse = []
    a_vect = []
    t_acc = []
    inertie = []
    for i in range(len(t)-1):
        v_vect.append((x_vect[i+1]-x_vect[i])/(t[i+1]-t[i]))
        t_vitesse.append((t[i]+t[i+1])/2)
    for i in range(len(t_vitesse)-1):
        a_vect.append((v_vect[i+1]-v_vect[i])/(t_vitesse[i+1]-t_vitesse[i]))
        t_acc.append((t_vitesse[i]+t_vitesse[i+1])/2)
    for i in range(len(t_acc)):
        inertie.append(a_vect[i]*masse_bille)
    #fig = plt.figure(1,figsize = (16, 12))
    #plt.plot(t_acc,a_vect)
    #plt.show()
    fig = plt.figure(1,figsize = (16, 12))
    plt.plot(t,f_vect,'.',markersize = 20)
    plt.plot(t_acc,inertie,'.',markersize = 10)
    plt.xlabel("temps (s)")
    plt.ylabel("Force (pN)")
    plt.show()
    return(t_relax) 

f0p = 40e-12
list_f0p = [15e-12,20e-12,25e-12,30e-12,35e-12,40e-12,45e-12,50e-12]
list_f0p_l_o_p = [log_ou_pas(15),log_ou_pas(20),log_ou_pas(25),log_ou_pas(30),log_ou_pas(35),log_ou_pas(40),log_ou_pas(45),log_ou_pas(50)]
print(len(list_f0p),len(list_f0p_l_o_p))
fig = plt.figure(1,figsize = (16, 12))
for coeff in [0.1,0.15,0.2,0.25,0.3,0.5]:
    list_t_relax = []
    for f0p in list_f0p:
        list_t_relax.append(resoudre_Fx_v(f0p,coeff))
    
    [absi,a,b] = fiit(list_f0p_l_o_p,list_t_relax)
    print(f'a = {a[0][0]}, b = {b[0]}')

    a = a[0][0]
    b = b[0]
    reg_force = list_f0p_l_o_p
    reg_time = [a*r_f+b for r_f in reg_force]
    plt.plot(reg_force,reg_time,'-.',color='grey')  
    plt.plot(list_f0p_l_o_p,list_t_relax,'.',label = str(coeff))
plt.legend()
plt.show()


eta = 10**(-3)
f0p = 40*10**(-9)
Lp = 56*10**(-9)
Lc = 16.5*10**(-6)
R = 4.98*10**(-6)
a = 2*10**(-9)
print(a)
pi = 3.141592653
xi_bille = 6*pi*eta*R

xi_ADN = 2*pi*eta/math.log(Lc/a,10)
print(xi_bille,xi_ADN)
t0 = xi_bille/4*(4.11/Lp)**(1/2)*Lc**2/f0p**(3/2)
print(t0)
  


# In[60]:


resoudre_Fx_v()


# In[ ]:





# In[ ]:




