import multiprocessing
import time
import main
import matplotlib.pyplot as plt

num_runs=10
num_cores=10

nc=min([num_runs,num_cores])

if __name__== "__main__":
    obliquity=2 # 1==0 deg; 2==45 deg; 3==90 deg
    sputtering=True
    photochemical_escape=True 
    pool = multiprocessing.Pool(processes=num_cores)
    results=[pool.apply_async(main.run_model, \
       args=([i],obliquity), \
       kwds={'sputtering_flag': sputtering,\
       'pce_flag':photochemical_escape}) for i in range(num_runs)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()
    

    fig=plt.figure()
    plt.ion()
    ax1=plt.subplot(211)
    COUNT=0
    for i in range(len(output)):
        (t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
        plt.plot(t, ystore[:,0])
        if (ystore[-1,0] < 0.1):
            COUNT=COUNT+1

    plt.yscale('log')
    plt.ylabel('Pressure (atm)')
    plt.text(0.5,0.5,'(a) Pressure: ' + str(COUNT) + ' collapsed',transform=ax1.transAxes)
    if obliquity == 1:
        plt.title('Obliquity = 0$^\circ$')
    elif obliquity == 2:
        plt.title('Obliquity = 45$^\circ$')
    elif obliquity == 3:
        plt.title('Obliquity = 90$^\circ$')

    ax2=plt.subplot(212)
    for i in range(len(output)):
    	(t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
    	plt.plot(t, ystore[:,1])
    plt.xlabel('time (Gy)')
    plt.ylabel(u'$\delta ^{15}$N [‰]')
    plt.text(0.5,0.2,'(b)',transform=ax2.transAxes)

    fig.tight_layout()
    
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

    

