import multiprocessing
import time
import main
import matplotlib.pyplot as plt
from isotopic_data import errorbars1
import numpy as np
num_runs=10
num_cores=10

nc=min([num_runs,num_cores])

if __name__== "__main__":
    
    obliquity=2 # 1==0 deg; 2==45 deg; 3==90 deg
    sputtering=True
    photochemical_escape=True 
    textadd=', with Sputtering and PCE'
    pool = multiprocessing.Pool(processes=num_cores)
    results=[pool.apply_async(main.run_model, \
       args=([i],obliquity), \
       kwds={'sputtering_flag': sputtering,\
       'pce_flag':photochemical_escape, \
       'C_Ne_IDP': 1., \
       'f_comet':0.0005, \
       'X_gas': 0.01, \
       'crater_model' : 1, \
       'dynamo_time' : 4.5}) for i in range(num_runs)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()
    

    fig=plt.figure(figsize=(8,15))
    plt.ion()
    ax1=plt.subplot(611)
    COUNT=0
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            COUNT=COUNT+1
            plt.plot(t, ystore[:,0])
        else:
            plt.plot(t, ystore[:,0],'-k',lw=0.1)

    plt.yscale('log')
    plt.ylabel('Pressure (atm)')
    plt.text(0.5,0.5,'(a) Pressure: ' + str(COUNT) + ' collapsed',transform=ax1.transAxes)
    if obliquity == 1:
        plt.title('Obliquity = 0$^\circ$' + textadd)
    elif obliquity == 2:
        plt.title('Obliquity = 45$^\circ$' + textadd)
    elif obliquity == 3:
        plt.title('Obliquity = 90$^\circ$' + textadd)
    elif obliquity == -1:
        plt.title('Control' + textadd)

    ax2=plt.subplot(612)
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            plt.plot(t, ystore[:,1])
        else:
            plt.plot(t, ystore[:,1],'-k',lw=0.1)
    #plt.xlabel('time (Gy)')
    plt.ylabel(u'$\delta ^{15}$N [‰]')
    plt.text(0.5,0.2,'(b)',transform=ax2.transAxes)
    ax2.errorbar([0],errorbars1[0][1],np.array([[errorbars1[0][1]-errorbars1[0][0],\
        errorbars1[0][2]-errorbars1[0][1]]]).T,np.array([[0,0]]).T,'r.',capsize=10) 

    ax3=plt.subplot(613)
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            plt.plot(t, 1.0/ystore[:,2])
        else:
            plt.plot(t, 1.0/ystore[:,2],'-k',lw=0.1)
    #plt.xlabel('time (Gy)')
    plt.ylabel(u'$\\frac{^{20}Ne}{^{22}Ne}$')
    plt.text(0.5,0.2,'(c)',transform=ax3.transAxes)
    ax3.errorbar([0],errorbars1[1][1],np.array([[errorbars1[1][1]-errorbars1[1][0],\
        errorbars1[1][2]-errorbars1[1][1]]]).T,np.array([[0,0]]).T,'r.',capsize=10) 

    ax4=plt.subplot(614)
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            plt.plot(t, ystore[:,3])
        else:
            plt.plot(t, ystore[:,3],'-k',lw=0.1)
    #plt.xlabel('time (Gy)')
    plt.ylabel(u'$\\frac{^{38}Ar}{^{36}Ar}$')
    plt.text(0.1,0.8,'(d)',transform=ax4.transAxes)
    ax4.errorbar([0],errorbars1[2][1],np.array([[errorbars1[2][1]-errorbars1[2][0],\
        errorbars1[2][2]-errorbars1[2][1]]]).T,np.array([[0,0]]).T,'r.',capsize=10) 

    ax5=plt.subplot(615)
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            plt.plot(t, ystore[:,4])
        else:
            plt.plot(t, ystore[:,4],'-k',lw=0.1)
    #plt.xlabel('time (Gy)')
    plt.ylabel(u'$\\frac{^{86}Kr}{^{84}Kr}$')
    plt.text(0.5,0.8,'(e)',transform=ax5.transAxes)
    ax5.errorbar([0],errorbars1[3][1],np.array([[errorbars1[3][1]-errorbars1[3][0],\
        errorbars1[3][2]-errorbars1[3][1]]]).T,np.array([[0,0]]).T,'r.',capsize=10) 

    ax6=plt.subplot(616)
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            plt.plot(t, ystore[:,5])
        else:
            plt.plot(t, ystore[:,5],'-k',lw=0.1)
    plt.xlabel('time (yr)')
    plt.ylabel(u'$\\frac{^{136}Xe}{^{130}Xe}$')
    plt.text(0.1,0.2,'(f)',transform=ax6.transAxes)
    ax6.errorbar([0],errorbars1[4][1],np.array([[errorbars1[4][1]-errorbars1[4][0],\
        errorbars1[4][2]-errorbars1[4][1]]]).T,np.array([[0,0]]).T,'r.',capsize=10) 

    #fig.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('/tmp/temp1.png')


  
    # mole_elements
    fig=plt.figure()
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        if (ystore[-1,0] < 0.1):
            COUNT=COUNT+1
            plt.plot(main.isotopic_data.element_mass,  \
				mole_elements/main.isotopic_data.solar_abundances/1e16,'-s',lw=1)
        else:            
            plt.plot(main.isotopic_data.element_mass,  \
				mole_elements/main.isotopic_data.solar_abundances/1e16,'-s',lw=0.1)
			
    plt.plot(main.isotopic_data.element_mass,  \
		main.isotopic_data.abundances1['Obs'][1],'r-s',linewidth=3)

    plt.yscale('log')
    plt.grid()
    plt.ylabel('Abundances wrt solar values')
    plt.xlabel('Mass number')
    plt.title('Elemental Abundances ' + textadd)
    plt.savefig('/tmp/temp2.png')

    

