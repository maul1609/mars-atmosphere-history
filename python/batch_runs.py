import multiprocessing
import time
import main

num_runs=5
num_cores=5

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
    pool = multiprocessing.Pool(processes=num_cores)
    results=[pool.apply_async(main.run_model, args=([i])) for i in range(num_runs)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()
