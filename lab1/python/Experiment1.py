import numpy as np
import matplotlib.pyplot as plt
import subprocess
from datetime import datetime
import argparse

def ImportDataTimeFileContent(file_name: str) -> list:
    res_list = []
    fd = open(file_name, "r")
    # the last symbol is \n
    # so, the last split str will be ""
    line = fd.readline()
    res = float(line)
    fd.close()
    return res

def CleanFile(file_path):
    fd = open(file_path, "w")
    fd.write("")
    fd.close()

def main():
    points_list = [100,   500,   1000,  2000, 4000, 6000, 8000, 12000, 16000, 20000]

    consistent_prog_time = []
    parallel_prog_time = []

    print("Start 1 - thread")
    for points in points_list:
        print(f"{points} points")
        subprocess.run(["../non_par_lab", f'{points}', f'{points}'])
        CleanFile("./data.txt")
            
        res = ImportDataTimeFileContent("time.txt")
        consistent_prog_time.append(np.mean(res))
        CleanFile("./time.txt")
    print("End 1 - thread")

    for i in range(2, 5):
        print(f"Start {i} - thread")
        i_thread_res = []
        for points in points_list:
            print(f"{points} points")
            subprocess.run(["mpiexec", "-n", f'{i}', "../lab", f'{points}', f'{points}'])
            CleanFile("./data.txt")

            res = ImportDataTimeFileContent("time.txt")
            i_thread_res.append(res)
            CleanFile("./time.txt")
        parallel_prog_time.append(i_thread_res)
        print(f"End {i} thread")
        

    date = datetime.strftime(datetime.now(), "%d.%m.%Y-%H.%M.%S")
    save_file_name = r"../images/" + date + r".jpg"

    fig = plt.figure()
    plt.scatter(np.array(points_list), np.array(consistent_prog_time))
    plt.plot(np.array(points_list), np.array(consistent_prog_time), label='1 - поток')
    plt.scatter(np.array(points_list), np.array(parallel_prog_time[0]))
    plt.plot(np.array(points_list), np.array(parallel_prog_time[0]), label='2 - поток')
    plt.scatter(np.array(points_list), np.array(parallel_prog_time[1]))
    plt.plot(np.array(points_list), np.array(parallel_prog_time[1]), label='3 - поток')
    plt.scatter(np.array(points_list), np.array(parallel_prog_time[1]))
    plt.plot(np.array(points_list), np.array(parallel_prog_time[2]), label='4 - поток')
    plt.xlabel("Количество точек")
    plt.ylabel("Время работы (мс)")
    plt.legend()
    plt.grid()
    fig.savefig(save_file_name)

if __name__ == '__main__':
    main()