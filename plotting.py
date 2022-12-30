import matplotlib.pyplot as plt
import numpy as np


cond_length=np.arange(0.075,0.201,0.001)    #0.075,0.2,0.001
cond_width=np.arange(0.08,0.17,0.01)
I0 = []#output current
I1 = []#cond current

z=[[1]* len(cond_width) for i in range(len(cond_length))]

results = open("results.txt","r")
for i in results:
    buffer = i.split(" ")
    I0.append(float(buffer[-2]))
    I1.append(float(buffer[-1]))

for i in range(len(cond_length)):
    for j in range(len(cond_width)):
        print(i)
        z[i][j] = I0[j+i*8]
print(z)
# fig5 = plt.figure()
# h5 = plt.contourf(cond_width,cond_length,z,10,cmap='jet')
# plt.title("output current depending on parameters")
# cbar = fig5.colorbar(h5,fraction=0.01)
# plt.xlabel("Capacitor width, mm")
# plt.ylabel("Capacitor length, mm")
# plt.show()


