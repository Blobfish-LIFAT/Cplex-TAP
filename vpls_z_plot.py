import os
import subprocess
from matplotlib import pyplot as plt

log = "C:\\Users\\chanson\\Desktop\\log_12_500_120s_8w.txt"
optimal = 99.3478

#log = "C:\\Users\\chanson\\Desktop\\log_22_500_120s_8w.txt"
#optimal = 98.978

clk_start = 0
clk_rate = 0
values = []
absolutes = []
times = []
iterations = []
vpls_started = False


with open(log) as file:
    for line in file:
        line = line.strip()
        if "CLK_RATE" in line:
            clk_rate = float(line.split(' ')[1])
        if "CLK_START" in line and "CLK_START_ITER" not in line:
            clk_start = float(line.split(' ')[1])
            vpls_started = True
        if "CLK_START_ITER" in line:
            iterations.append(float(line.split(' ')[1]))
        if "[INFO MIP Callback]" in line and vpls_started:
            tmp = line.split(" ")
            times.append((float(tmp[4]) - clk_start)/clk_rate)
            #values.append((float(tmp[6])/optimal)*100)
            values.append( ((optimal - float(tmp[6]))/optimal)*100 )
            absolutes.append(float(tmp[6]))



plt.scatter(times, values)
#plt.plot([], [], ' ', label="Opt Z =" + str(optimal))
#plt.axhline(y=optimal, label="Optimal " + str(optimal), ls='dotted')

for i in iterations:
    t = (i - clk_start)/clk_rate
    plt.axvline(x=t, ls='dashed', color="red")


plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Objective value (% from optimal)", fontsize=14)
plt.legend()
plt.show()

print("iteration,time (absolute),time (relative),gap to optimal")
for i, tstart in enumerate(iterations):
    tstart = (tstart - clk_start)/clk_rate
    for j, time in enumerate(times):
        if time > tstart:
            if len(iterations) - 1 > i and time > (iterations[i+1]- clk_start)/clk_rate:
                continue
            if absolutes[j] > absolutes[j-1]:
                print(",".join(map(str, [i, time, time - tstart, values[j]])))