import multiprocessing
import time
import main
import matplotlib.pyplot as plt

num_runs=20
num_cores=10

nc=min([num_runs,num_cores])

if __name__== "__main__":
#    manager = multiprocessing.Manager()
#    return_dict = manager.dict()
# 	jobs = []
# 	for i in range(num_runs):
# 		p = multiprocessing.Process(target=main.run_model, args=(i, return_dict))
# 		jobs.append(p)
# 		p.start()
# 	
# 	for proc in jobs:
# 		proc.join()
# 		
    obliquity=1 # 1==0 deg; 2==45 deg; 3==90 deg
    pool = multiprocessing.Pool(processes=num_cores)
    results=[pool.apply_async(main.run_model, args=([i],obliquity)) for i in range(num_runs)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()
    

    plt.figure()
    plt.ion()
    plt.subplot(211)
    for i in range(len(output)):
    	(t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
    	plt.plot(t, ystore[:,0])
    	plt.yscale('log')
    plt.subplot(212)
    for i in range(len(output)):
    	(t,ystore,mole_elements,isotopes_sim,Amars)=output[i]
    	plt.plot(t, ystore[:,1])