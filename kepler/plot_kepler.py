import pyart
import cartopy.crs as ccrs

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.patheffects as pe

# from matplotlib.colors import ListedColormap
import json
import sys
from utils import wescon_kml_grid
import os
import datetime as dt


variables = ["DBZ", "LDR", "VEL"]
desired_elevations = [1.0, 1.5, 2.5]


def make_ppi_plots(
    radar,
    radar_name,
    outdir,
    var,
    var_scales,
    outer_lines,
    h_lines,
    v_lines,
    storm_boxes,
    labels,
    xlim=[-4.6, -1],
    ylim=[49.75, 52.25],
    desired_elevations=[1.0],
):
    # get variable info
    vmin = var_scales[radar_name][var]["min"]
    vmax = var_scales[radar_name][var]["max"]
    # num_colours = var_scales[radar_name][var]["num_colours"]
    colourmap = var_scales[radar_name][var]["colourmap"]

    # radar display
    display = pyart.graph.RadarMapDisplay(radar)

    # make plots for each sweep
    for sweep in range(radar.nsweeps):
        sweep_elevation = radar.elevation["data"][
            sum(radar.rays_per_sweep["data"][0 : sweep + 1]) - 1
        ]
        if float(sweep_elevation) in desired_elevations:
            if not os.path.exists(f"{outdir}/{sweep_elevation}deg"):
                os.mkdir(f"{outdir}/{sweep_elevation}deg")
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
            display.plot_ppi_map(
                var,
                ax=ax,
                resolution="10m",
                vmin=vmin,
                vmax=vmax,
                sweep=sweep,
                colorbar_orient="horizontal",
                cmap=colourmap,
                mask_tuple=("SNR",1),
            )
            display.plot_range_rings(
                [10,20,30,40],
                col='black',
                lw=1,
            )
            gl = ax.gridlines(
                crs=ccrs.PlateCarree(),
                draw_labels=True,
                linewidth=1,
                color="gray",
                alpha=0.3,
                linestyle="--",
            )
            gl.top_labels = False
            gl.right_labels = False
            """ 
            for line in outer_lines:
                xs = [
                    float(line[0].split(" ")[0].split(",")[0]),
                    float(line[0].split(" ")[1].split(",")[0]),
                ]
                ys = [
                    float(line[0].split(" ")[0].split(",")[1]),
                    float(line[0].split(" ")[1].split(",")[1]),
                ]
                ax.add_artist(lines.Line2D(xs, ys, color="black"))
            for line in h_lines:
                xs = [
                    float(line[0].split(" ")[0].split(",")[0]),
                    float(line[0].split(" ")[1].split(",")[0]),
                ]
                ys = [
                    float(line[0].split(" ")[0].split(",")[1]),
                    float(line[0].split(" ")[1].split(",")[1]),
                ]
                ax.add_artist(lines.Line2D(xs, ys, color="black"))
            for line in v_lines:
                xs = [
                    float(line[0].split(" ")[0].split(",")[0]),
                    float(line[0].split(" ")[1].split(",")[0]),
                ]
                ys = [
                    float(line[0].split(" ")[0].split(",")[1]),
                    float(line[0].split(" ")[1].split(",")[1]),
                ]
                ax.add_artist(lines.Line2D(xs, ys, color="black"))
            for line in storm_boxes:
                xs = [float(line[i].split(",")[0]) for i in range(len(line))]
                ys = [float(line[i].split(",")[1]) for i in range(len(line))]
                xy = (min(xs), min(ys))
                width = max(xs) - min(xs)
                height = max(ys) - min(ys)
                ax.add_artist(
                    mpl.patches.Rectangle(
                        xy,
                        width,
                        height,
                        facecolor="#FFFFFF00",
                        edgecolor="blue",
                        linewidth=2,
                    )
                )
            for k, v in labels.items():
                x = float(str(v).split(",")[0])
                y = float(str(v).split(",")[1])
                txt = mpl.text.Text(
                    x=x,
                    y=y,
                    text=k,
                    horizontalalignment="center",
                    verticalalignment="center",
                )
                txt.set_path_effects([pe.withStroke(linewidth=2, foreground="w")])
                ax.add_artist(txt)
            """
            plt.xlim(xlim)
            plt.ylim(ylim)

            sweep_time = dt.datetime.strftime(
                pyart.graph.common.generate_radar_time_sweep(radar, sweep),
                "%Y%m%d%H%M%S",
            )
            ax.set_aspect(1./ax.get_data_ratio())
            plt.savefig(f"{outdir}/{sweep_elevation}deg/{sweep_time}00{var}.ppi.png")
            plt.close()


def make_rhi_plots(radar, radar_name, outdir, var, var_scales, ylim=[0, 12]):
    # make out dir
    if not os.path.exists(f"{outdir}/rhi"):
        os.mkdir(f"{outdir}/rhi")
    # get variable info
    vmin = var_scales[radar_name][var]["min"]
    vmax = var_scales[radar_name][var]["max"]
    # num_colours = var_scales[radar_name][var]["num_colours"]
    colourmap = var_scales[radar_name][var]["colourmap"]

    # radar display
    display = pyart.graph.RadarMapDisplay(radar)

    # make plots
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    display.plot_rhi(
        var,
        ax=ax,
        vmin=vmin,
        vmax=vmax,
        colorbar_orient="horizontal",
        cmap=colourmap,
        mask_tuple=("SNR",1),
    )
    plt.grid(linewidth=1, color="gray", alpha=0.3, linestyle="--")
    plt.ylim(ylim)
    sweep_time = dt.datetime.strftime(
        pyart.graph.common.generate_radar_time_sweep(radar, 0), "%Y%m%d%H%M%S"
    )
    plt.savefig(f"{outdir}/rhi/{sweep_time}00{var}.rhi.png")
    plt.close()


def main(
    radarfile, outdir, variables=variables, radar_name="ncas-radar-mobile-ka-band-1"
):
    file_path = "/".join(__file__.split("/")[:-1])
    radar = pyart.io.read(radarfile)

    # read in stuff from utils
    outer_lines, h_lines, v_lines, storm_boxes, labels = wescon_kml_grid.read_kml(
        f"{file_path}/../WesConGrid.kml"
    )
    with open(f"{file_path}/utils/var_scales.json") as f:
        var_scales = json.load(f)

    # what type of radar scan
    if "ppi" in radarfile:
        scan_type = "PPI"
    elif "rhi" in radarfile:
        scan_type = "RHI"
    else:
        print("Unknown scan type, quitting...")
        sys.exit()

    for var in variables:
        # make plot
        if scan_type == "PPI":
            make_ppi_plots(
                radar,
                radar_name,
                outdir,
                var,
                var_scales,
                outer_lines,
                h_lines,
                v_lines,
                storm_boxes,
                labels,
                desired_elevations=desired_elevations,
                xlim=[-2.6,-1.4],
                ylim=[51.1,51.9],
            )
        elif scan_type == "RHI":
            make_rhi_plots(
                radar,
                radar_name,
                outdir,
                var,
                var_scales,
            )


if __name__ == "__main__":
    radarfile = sys.argv[1]
    outdir = sys.argv[2]
    main(radarfile, outdir)
