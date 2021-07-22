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
    taille = 10

    #plt.scatter(times, values, s=taille, color="black", label="CPLEX")

    log = "logs/ist12_best_procedural.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=taille, color="green", label="Left/Right")

    log = "logs/ist12_best_random.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=taille, color="red", label="Random")

    log = "logs/hamming_12.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=taille, color="orange", label="Hamming S")

    log = "logs/hamming_sx_ist12.log"
    iterations, values, clk_start, clk_rate, times, absolutes = read_log(log, optimal)
    plt.scatter(times, values, s=taille, color="blue", label="Hamming SX")

    plt.xlabel("Time (s)", fontsize=14)
    plt.ylabel("Objective value (% from optimal)", fontsize=14)
    plt.ylim(ymin=0)
    
    plt.title("Istance 12")

    major_ticks = np.arange(0, np.max(times)+1, 100)
    minor_ticks = np.arange(0, np.max(times)+1, 50)
    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(np.arange(0, np.max(values) + 1, 1), minor=True)

    plt.legend(fontsize=18)
    plt.show()

