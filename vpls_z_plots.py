import os
import subprocess
import numpy as np
from matplotlib import pyplot as plt

def read_log(log, optimal):
    clk_start = 0
    clk_rate = 0
    values = []
    absolutes = []
    times = []
    iterations = []
    vpls_started = False
    params = ""
    cplex_init = 0
    with open(log) as file:
        lnum = 0
        for line in file:
            line = line.strip()
            if "PARAMS:" in line:
                tmp = line.split(" ")
                params = "Ep. " + tmp[1] + "/" + tmp[2] + " | h=" + tmp[3] + " | times " + tmp[4] + "/" + tmp[5]
                cplex_init = int(tmp[4])
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
                values.append( ((optimal - float(tmp[6]))/optimal)*100 )
                absolutes.append(float(tmp[6]))
            lnum += 1
    return iterations, values, clk_start, clk_rate, times, absolutes

if __name__ == '__main__':

    log = "logs/ist12_full.log"
    optimal = 99.3478
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    plt.scatter(times, values, s=5, color="black", label="CPLEX")

    log = "logs/procedural_45_8.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=5, color="green", label="8")

    log = "logs/procedural_45_10.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=5, color="red", label="10")

    log = "logs/procedural_45_12.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=5, color="orange", label="12")

    log = "logs/procedural_45_14.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=5, color="blue", label="14")

    log = "logs/procedural_45_16.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=5, color="purple", label="16")

    plt.xlabel("Time (s)", fontsize=14)
    plt.ylabel("Objective value (% from optimal)", fontsize=14)
    plt.ylim(ymin=0)
    
    plt.title("Ist 12 - 45s")

    major_ticks = np.arange(0, np.max(times)+1, 100)
    minor_ticks = np.arange(0, np.max(times)+1, 50)
    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(np.arange(0, np.max(values) + 1, 1), minor=True)

    plt.legend()
    plt.show()

