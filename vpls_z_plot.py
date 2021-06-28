import os
import subprocess
import numpy as np
from matplotlib import pyplot as plt

log = "logs/hamming_12.log"

#OPT 12/500
optimal = 99.3478


clk_start = 0
clk_rate = 0
values = []
absolutes = []
times = []
iterations = []
vpls_started = False
params = ""
conv = "max"
stop_type = []
solutions = []
cplex_init = 0

with open(log) as file:
    lnum = 0
    for line in file:
        line = line.strip()
        if "PARAMS:" in line:
            tmp = line.split(" ")
            params = "Ep. " + tmp[1] + "/" + tmp[2] + " | h=" + tmp[3] + " | times " + tmp[4] + "/" + tmp[5]
            cplex_init = int(tmp[4])
        if vpls_started and "SOLUTION: " in line:
            solutions.append(line.replace("SOLUTION: ", ""))
        if vpls_started and "Status:" in line:
            stop_type.append(line.split("Status:")[1].strip())
            print(lnum)
        if vpls_started and "VPLS Converged at iteration" in line:
            conv = line.split("VPLS Converged at iteration")[1]
        if "CLK_RATE" in line:
            clk_rate = float(line.split(' ')[1])
        if "CLK_START" in line and "CLK_START_ITER" not in line:
            clk_start = float(line.split(' ')[1])
            vpls_started = True
        if "CLK_START_ITER" in line:
            iterations.append(float(line.split(' ')[1]))
        if vpls_started and "[INFO MIP Callback]" in line:
            tmp = line.split(" ")
            times.append((float(tmp[4]) - clk_start)/clk_rate)
            #values.append((float(tmp[6])/optimal)*100)
            values.append( ((optimal - float(tmp[6]))/optimal)*100 )
            absolutes.append(float(tmp[6]))
        lnum += 1


print("iteration,time (absolute),time (relative),gap to optimal")
for i, tstart in enumerate(iterations):
    tstart = (tstart - clk_start)/clk_rate
    for j, time in enumerate(times):
        if time > tstart:
            if len(iterations) - 1 > i and time > (iterations[i+1]- clk_start)/clk_rate:
                continue
            if absolutes[j] > absolutes[j-1]:
                print(",".join(map(str, [i, time, time - tstart, values[j]])))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

#plt.plot([], [], ' ', label="Opt Z =" + str(optimal))
#plt.axhline(y=optimal, label="Optimal " + str(optimal), ls='dotted')

for c, i in enumerate(iterations[1:]):
    if c == 0:
        continue
    t = (i - clk_start)/clk_rate
    color = ''
    if stop_type[c] == "Feasible":
        color = "blue"
    elif stop_type[c] == "Optimal":
        color = "green"
    if c > 0 and solutions[c-1] == solutions[c] and color == "green":
        color = "black"
    elif c > 0 and solutions[c-1] == solutions[c] and color == "blue":
        color = "red"
    plt.axvline(x=t + cplex_init, ls='dashed', color=color)

plt.scatter(times, values, s=5, color="black")

plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Objective value (% from optimal)", fontsize=14)
plt.title(params + " | Conv. " + conv + " | Zf gap=" + str(values[-1]))
plt.ylim(ymin=0)

major_ticks = np.arange(0, np.max(times)+1, 100)
minor_ticks = np.arange(0, np.max(times)+1, 50)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(np.arange(0, np.max(values) + 1, 1), minor=True)

plt.legend()
plt.show()

