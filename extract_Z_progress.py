import os
import subprocess
from matplotlib import pyplot as plt


args = ("logs/res3.log")
popen = subprocess.Popen(args, stdout=subprocess.PIPE)

clk_start = 0
clk_rate = 0
values = []
times = []
optimal = 0
total = "None"

for line in popen.stdout:
    line = str(line, encoding="utf-8").strip()
    if "CLK_RATE" in line:
        clk_rate = float(line.split(' ')[1])
    if "CLK_START" in line:
        clk_start = float(line.split(' ')[1])
    if "[INFO MIP Callback]" in line:
        tmp = line.split(" ")
        times.append((float(tmp[4]) - clk_start)/clk_rate)
        values.append(float(tmp[6]))
        #print(times[-1], values[-1])
    if "Objective:" in line:
        optimal = float(line.split('ve: ')[1])
    if "TIME TO SOLVE " in line:
        total = line


plt.scatter(times, values)

plt.plot([], [], ' ', label="Opt Z =" + str(optimal))
plt.plot([], [], ' ', label=total)

for i, v in enumerate(values):
    if v > optimal - optimal*0.01:
        plt.axvline(x=times[i], label="Time to 1% of optimal " + str(times[i]) + " s", ls='--')
        break

for i, v in enumerate(values):
    if v == optimal:
        plt.axvline(x=times[i], label="Time to optimal " + str(times[i]) + " s", ls='--')
        break

plt.xlabel("Time (s)", fontsize=14)
plt.ylabel("Objective value", fontsize=14)
plt.legend()

plt.savefig("/tmp/out_plt.png")