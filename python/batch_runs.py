import multiprocessing
import time
import main

num_runs=100
num_cores=100

nc=min([num_runs,num_cores])

if __name__== "__main__":
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
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
    for i in range(num_runs):
        pool.apply_async(main.run_model, args=(i,return_dict))
    pool.close()
    pool.join()
