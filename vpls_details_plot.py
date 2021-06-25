import os
import subprocess
import numpy as np
from matplotlib import pyplot as plt

log = "logs/procedural_1.log"

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
            params = "Ep. " + tmp[1] + "/" + tmp[2] + " | h=" + tmp[3]
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

selected = 10

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

sel_times = []
sel_values = []
for i, tstart in enumerate(iterations):
    tstart = (tstart - clk_start)/clk_rate
    if i == selected:
        for j, time in enumerate(times):
            if time > tstart:
                if len(iterations) - 1 > i and time > (iterations[i+1]- clk_start)/clk_rate:
                    continue
                sel_times.append(time - tstart)
                sel_values.append(values[j])


plt.scatter(sel_times, sel_values, s=5, color="black")

plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Objective value (% from optimal)", fontsize=14)
plt.title(params + " | Iteration = " + str(selected))
plt.ylim(ymin=0)

plt.legend()
plt.show()

