import re
import numpy as np
import matplotlib.pyplot as plt
import json
import sys


def read_ascii_file(filename):
    with open(filename, encoding="charmap") as f:
        data = f.readlines()

    data = [i.strip() for i in data]
    headers = None
    time = None
    i = 0

    while not headers and i < len(data):
        if re.match("\d{2}/\d{2}/\d{2}\s\d{2}:\d{2}", data[i]):
            file_time = data[i].split(" ")[1]
            hours, minutes = file_time.split(":")
            time = (int(hours) * 60) + int(minutes)
            i += 1
        elif data[i].startswith("Altitude") and data[i].endswith("SNR"):
            headers = data[i]
        else:
            i += 1

    while "  " in headers:
        headers = headers.replace("  ", " ")
    headers = headers.split(" ")

    if not headers == ["Altitude", "East", "North", "vertical", "dd", "ff", "SNR"]:
        print("Trouble")
        sys.exit()

    altitudes = []
    east_wind = []
    north_wind = []
    vertical_wind = []
    wind_dir = []
    wind_speed = []
    snr = []

    i = i + 2  # start of data
    data_end = False

    while not data_end and i < len(data):
        if len(data[i]) > 0:
            line = data[i]
            while "  " in line:
                line = line.replace("  ", " ")
            line = line.split(" ")
            if line[0] == "*":
                line = line[1:]
            line = [np.nan if d == "//////" else float(d) for d in line]
            altitudes.append(line[0])
            east_wind.append(line[1])
            north_wind.append(line[2])
            vertical_wind.append(line[3])
            wind_dir.append(line[4])
            wind_speed.append(line[5])
            snr.append(line[6])
            i += 1
        else:
            data_end = True

    return (
        time,
        altitudes,
        east_wind,
        north_wind,
        vertical_wind,
        wind_dir,
        wind_speed,
        snr,
    )


def read_all_files(infiles):
    time_arr = np.ones(len(infiles)) * -1
    altitudes_arr = np.ones(len(infiles)) * -1
    east_wind_arr = np.ones(len(infiles)) * -1
    north_wind_arr = np.ones(len(infiles)) * -1
    vertical_wind_arr = np.ones(len(infiles)) * -1
    wind_dir_arr = np.ones(len(infiles)) * -1
    wind_speed_arr = np.ones(len(infiles)) * -1
    snr_arr = np.ones(len(infiles)) * -1
    for i, infile in enumerate(infiles):
        (
            time,
            altitudes,
            east_wind,
            north_wind,
            vertical_wind,
            wind_dir,
            wind_speed,
            snr,
        ) = read_ascii_file(infile)
        if i == 0:
            altitudes_arr = np.ones((len(infiles), len(altitudes))) * -1
            east_wind_arr = np.ones((len(infiles), len(altitudes))) * -1
            north_wind_arr = np.ones((len(infiles), len(altitudes))) * -1
            vertical_wind_arr = np.ones((len(infiles), len(altitudes))) * -1
            wind_dir_arr = np.ones((len(infiles), len(altitudes))) * -1
            wind_speed_arr = np.ones((len(infiles), len(altitudes))) * -1
            snr_arr = np.ones((len(infiles), len(altitudes))) * -1
        time_arr[i] = time
        altitudes_arr[i] = altitudes
        east_wind_arr[i] = east_wind
        north_wind_arr[i] = north_wind
        vertical_wind_arr[i] = vertical_wind
        wind_dir_arr[i] = wind_dir
        wind_speed_arr[i] = wind_speed
        snr_arr[i] = snr
    return (
        time_arr,
        altitudes_arr,
        east_wind_arr,
        north_wind_arr,
        vertical_wind_arr,
        wind_dir_arr,
        wind_speed_arr,
        snr_arr,
    )


def plot_variable(
    var,
    time,
    altitude,
    outdir=".",
    filename="plot",
    vmin=None,
    vmax=None,
    cmap="viridis",
):
    time_plot = time
    for i in range(np.shape(altitude)[1] - 1):
        time_plot = np.vstack([time, time_plot])
    c = plt.pcolormesh(time_plot.T, altitude, var, vmin=vmin, vmax=vmax, cmap=cmap)
    plt.colorbar(c)
    plt.savefig(f"{outdir}/{filename}.png")
    plt.close()


def plot_winds(
    east_wind,
    north_wind,
    time,
    altitude,
    outdir=".",
    filename="wind_plot",
    vmin=None,
    vmax=None,
    cmap="viridis",
    barb_interval=4,
):
    time_plot = time
    for i in range(np.shape(altitude)[1] - 1):
        time_plot = np.vstack([time, time_plot])

    wind_speed = ((east_wind**2) + (north_wind**2)) ** 0.5
    # want all arrows to be same length
    arrow_length = 2
    v1 = (
        ((arrow_length**2) / ((east_wind**2 / north_wind**2) + 1)) ** 0.5
    ) * np.sign(north_wind)
    u1 = (v1 / north_wind) * east_wind

    c = plt.pcolormesh(
        time_plot.T, altitude, wind_speed, vmin=vmin, vmax=vmax, cmap=cmap
    )
    plt.colorbar(c)

    plt.quiver(
        time_plot.T[::barb_interval, ::barb_interval],
        altitude[::barb_interval, ::barb_interval],
        u1[::barb_interval, ::barb_interval],
        v1[::barb_interval, ::barb_interval],
        scale=48,
        scale_units="width",
    )
    plt.savefig(f"{outdir}/{filename}.png")
    plt.close()


def main(infiles, outdir):
    file_path = "/".join(__file__.split("/")[:-1])
    if file_path == "":
        file_path = "."
    with open(f"{file_path}/utils/var_scales.json") as f:
        var_scales = json.load(f)["ncas-radar-wind-profiler-2"]

    (
        time_arr,
        altitudes_arr,
        east_wind_arr,
        north_wind_arr,
        vertical_wind_arr,
        wind_dir_arr,
        wind_speed_arr,
        snr_arr,
    ) = read_all_files(infiles)

    plot_variable(
        snr_arr,
        time_arr,
        altitudes_arr,
        outdir=outdir,
        filename="nrwp2_snr",
        vmin=var_scales["SNR"]["min"],
        vmax=var_scales["SNR"]["max"],
        cmap=var_scales["SNR"]["colourmap"],
    )

    plot_variable(
        vertical_wind_arr,
        time_arr,
        altitudes_arr,
        outdir=outdir,
        filename="nrwp2_upward_wind",
        vmin=var_scales["upward_wind"]["min"],
        vmax=var_scales["upward_wind"]["max"],
        cmap=var_scales["upward_wind"]["colourmap"],
    )

    plot_winds(
        east_wind_arr,
        north_wind_arr,
        time_arr,
        altitudes_arr,
        outdir=outdir,
        filename="nrwp2_winds",
        vmin=var_scales["wind_speed"]["min"],
        vmax=var_scales["wind_speed"]["max"],
        cmap=var_scales["wind_speed"]["colourmap"],
    )


if __name__ == "__main__":
    outdir = sys.argv[1]
    infiles = sys.argv[2:]
    main(infiles, outdir)
