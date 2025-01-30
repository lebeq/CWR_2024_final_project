import numpy as np
from matplotlib import pyplot as plt

################################# Plot Set-Up f. Teil 4 ####################################

def plot_setup(file_name: str, N: int):
    S_0, mc_callprice, var, ana_callprice = np.loadtxt(file_name, float, delimiter=',',unpack=True)
    plt.errorbar(S_0,mc_callprice,yerr=var, color= 'red', fmt='o',label='Monte-Carlo')  
    plt.xlabel('Aktienwert zu Beginn ($S_0$) [EUR]')
    plt.ylabel('Preis [EUR]')
    plt.plot(S_0, ana_callprice, label='Analytische Lösung')
    plt.title(f'N = {N}')
    plt.grid()
    plt.minorticks_on()
    plt.legend()
    plt.savefig(f'N{N}_data.jpg')
    plt.show()

################################# Teil 1 ####################################################

##### ein Bsp für einen Aktienkursverlauf
t, S, ana_call_price = np.loadtxt('aktienkurs_data.csv', float, delimiter=',', unpack = True)
plt.plot(t,S)
plt.xlabel('Zeit [Jahre]')
plt.ylabel('Aktienkurs [EUR]')
plt.suptitle('Möglicher Verlauf eines Aktienkurses')
plt.grid()
plt.minorticks_on()
plt.savefig('./Aktienkurs.jpg')
plt.show()

####### Analytische Callprice Fkt i.A. von Anfangspreis der Aktie (S_0)
S_mat, new_ana_cp = np.loadtxt("analytic_solution.csv", float, delimiter=',',unpack=True)
plt.plot(S_mat, new_ana_cp, color='green')
plt.suptitle('Europäische Call-Option nach Black-Scholes')
plt.xlabel('Aktienwert zu Beginn ($S_0$) [EUR]')
plt.ylabel('Preis [EUR]')
plt.minorticks_on()
plt.grid()
plt.savefig('Analytischer_Callprice.jpg')
plt.show()

############################## Teil 2 #####################################################
x,D = np.loadtxt("aktien_pdf.csv", delimiter=',', unpack = True)
plt.bar(x, D, align='edge', width=(x[1]-x[0]), edgecolor='k')
plt.xlabel('S_T aus M-C')
plt.ylabel('Relative frequenz')
plt.title('PDF der Simulationsergebnisse')
plt.savefig('Teil3.png')
plt.show()

############################## Teil 4 Plots separat #######################################################

######### N = 100 #########
plot_setup("N100_data.csv",100)

######### N = 500 #########
plot_setup("N500_data.csv",500)

######### N = 1000 #########
plot_setup("N1000_data.csv",1000)

######### N = 10 000 #########
plot_setup("N10000_data.csv",10000)

########################## Teil 4 alle plots zusammen ##############################
fig, ax = plt.subplots(2,2, figsize = (10,8))
S_0_100, mc_callprice_100, var_100, ana_callprice_100 = np.loadtxt("N100_data.csv", float, delimiter=',',unpack=True)
S_0_500, mc_callprice_500, var_500, ana_callprice_500 = np.loadtxt("N500_data.csv", float, delimiter=',',unpack=True)
S_0_1000, mc_callprice_1000, var_1000, ana_callprice_1000 = np.loadtxt("N1000_data.csv", float, delimiter=',',unpack=True)
S_0_10000, mc_callprice_10000, var_10000, ana_callprice_10000 = np.loadtxt("N10000_data.csv", float, delimiter=',',unpack=True)

############ N = 100 ############
ax[0,0].errorbar(S_0_100,mc_callprice_100,yerr=var_100, color= 'red', fmt='o',label='Monte-Carlo')  
ax[0,0].plot(S_0_100, ana_callprice_100, label='Analytische Lösung')
ax[0,0].set_title('N = 100')
ax[0,0].grid()
ax[0,0].minorticks_on()
ax[0,0].legend()

########### N = 500 #############
ax[0,1].errorbar(S_0_500,mc_callprice_500,yerr=var_500, color= 'red', fmt='o',label='Monte-Carlo')  
ax[0,1].plot(S_0_500, ana_callprice_500, label='Analytische Lösung')
ax[0,1].set_title('N = 500')
ax[0,1].grid()
ax[0,1].minorticks_on()
ax[0,1].legend()

########## N = 1000 ###########
ax[1,0].errorbar(S_0_1000,mc_callprice_1000,yerr=var_1000, color= 'red', fmt='o',label='Monte-Carlo')  
ax[1,0].plot(S_0_1000, ana_callprice_1000, label='Analytische Lösung')
ax[1,0].set_title('N = 1000')
ax[1,0].grid()
ax[1,0].minorticks_on()
ax[1,0].legend()

########## N = 10 000 ##########
ax[1,1].errorbar(S_0_10000,mc_callprice_10000,yerr=var_10000, color= 'red', fmt='o',label='Monte-Carlo')  
ax[1,1].plot(S_0_10000, ana_callprice_10000, label='Analytische Lösung')
ax[1,1].set_title('N = 10 000')
ax[1,1].grid()
ax[1,1].minorticks_on()
ax[1,1].legend()

fig.supxlabel('Aktienwert zu Beginn ($S_0$) [EUR]')
fig.supylabel('Preis [EUR]')
fig.tight_layout()
plt.savefig('Teil4.jpg')
plt.show()

