import numpy as np
import scipy.interpolate as scint 

data= dict()

# figure 1, global annual mean surface temperature
# active clouds
data['active_clouds'] = np.array([[0.088328076,0.309148265,0.496845426,0.993690852,\
	1.965299685,3.003154574,3.831230284,4.659305994,6.403785489], \
    [199.4915254,208.220339,214.4915254,223.8983051,231.1016949,227.6271186,216.1864407,\
	218.8135593,225]])
# inactive clouds
data['inactive_clouds'] = np.array([[0.088328076,0.287066246,0.507886435,0.993690852,\
	1.976340694,2.914826498,3.798107256,4.648264984,6.348580442], \
    [199.0677966,204.9152542,209.5762712,215.4237288,219.9152542,219.8305085, \
    218.3050847,219.5762712,226.6949153]])
    

meanTvP1=scint.interp1d(data['active_clouds'][0],data['active_clouds'][1],fill_value='extrapolate',kind='linear')
meanTvP2=scint.interp1d(data['inactive_clouds'][0],data['inactive_clouds'][1],fill_value='extrapolate',kind='linear')

# figure 14, boundary for collapse (pressure vs obliquity)
data['collapse1'] = np.array([[0.8,0.893491442,1.004825929,1.15768754,1.266460627, \
	1.375916251,1.46923831,1.552717505,1.635281944,1.752237207,1.845410937, \
	1.950263512,1.997990345,2.068204665,2.10423179,2.111512177], \
    [29.94423792,29.73977695,29.22862454,28.00185874,26.9795539,25.65055762, \
	24.42379182,23.19702602,21.86802974,19.72118959,17.67657993,15.01858736,\
	13.38289963,10.31598513,6.942379182,5]])



data['collapse2']=4.452422411866442


if __name__ == "__main__":
	P=np.linspace(0,7,100)
	import matplotlib.pyplot as plt
	
	plt.ion()
	plt.plot(P,meanTvP1(P))
	plt.plot(P,meanTvP2(P))
	plt.show()
	