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
data['collapse1'] = np.array([[0.08,0.091865426,0.110467595,0.150790727,0.183426726,\
	0.223126214,0.304572251,0.392467945,0.425442056,0.511591388,0.574083387,0.594279167], \
    [30.01506024,29.68885853,29.53002733,28.54757949,27.560862,26.24281922,23.43808223,\
	20.13422246,18.81191001,14.34584989,9.048914731,5.0]])

meanPvObliquity=scint.interp1d(data['collapse1'][1],data['collapse1'][0],fill_value='extrapolate',kind='linear')


data['collapse2']=2.7521926650675255


if __name__ == "__main__":
	P=np.linspace(0,7,100)
	import matplotlib.pyplot as plt
	
	plt.ion()
	plt.plot(P,meanTvP1(P))
	plt.plot(P,meanTvP2(P))
	plt.xlabel('Mars surface pressure (bar)')
	plt.ylabel('Global mean surface temperature (K)')
	plt.legend(['Active clouds','Inactive clouds'])
	plt.show()
	