import pyart
import cartopy.crs as ccrs

# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as lines

# from matplotlib.colors import ListedColormap
import json
import sys
from utils import wescon_kml_grid
import os

variables = ["dBZ", "ZDR", "KDP", "RhoHV"]


def make_ppi_plots(
    radar,
    radar_name,
    outdir,
    var,
    var_scales,
    outer_lines,
    h_lines,
    v_lines,
    xlim=[-3.6, -1],
    ylim=[50.25, 51.75],
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
        if int(sweep_elevation) != 90:
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
            for line in outer_lines:
                xs = [float(line[0].split(",")[0]), float(line[1].split(",")[0])]
                ys = [float(line[0].split(",")[1]), float(line[1].split(",")[1])]
                ax.add_artist(lines.Line2D(xs, ys, color="black"))
            for line in h_lines:
                xs = [float(line[0].split(",")[0]), float(line[1].split(",")[0])]
                ys = [float(line[0].split(",")[1]), float(line[1].split(",")[1])]
                ax.add_artist(lines.Line2D(xs, ys, color="black"))
            for line in v_lines:
                xs = [float(line[0].split(",")[0]), float(line[1].split(",")[0])]
                ys = [float(line[0].split(",")[1]), float(line[1].split(",")[1])]
                ax.add_artist(lines.Line2D(xs, ys, color="black"))
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.savefig(
                f"{outdir}/{sweep_elevation}deg/{radar.metadata['start_datetime'].replace('-','').replace(':','').replace('Z','').replace('T','')}00{var}.ppi.png"
            )
            print(
                f"{outdir}/{sweep_elevation}deg/{radar.metadata['start_datetime'].replace('-','').replace(':','').replace('Z','').replace('T','')}00{var}.ppi.png written"
            )
            plt.close()


def make_rhi_plots(radar, radar_name, outdir, var, var_scales, ylim=[0, 15]):
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
    )
    plt.grid(linewidth=1, color="gray", alpha=0.3, linestyle="--")
    plt.ylim(ylim)
    plt.savefig(
        f"{outdir}/{radar_name}_{radar.metadata['start_datetime'].replace('-','').replace(':','')}_{var}_rhi.png"
    )
    print(
        f"{outdir}/{radar_name}_{radar.metadata['start_datetime'].replace('-','').replace(':','')}_{var}_rhi.png written"
    )
    plt.close()


def main(
    radarfile, outdir, variables=variables, radar_name="ncas-mobile-x-band-radar-2"
):
    file_path="/".join(__file__.split("/")[:-1])
    radar = pyart.io.read(radarfile)

    # read in stuff from utils
    outer_lines, h_lines, v_lines = wescon_kml_grid.read_kml(f"{file_path}/../WesConGrid.kml")
    with open(f"{file_path}/utils/var_scales.json") as f:
        var_scales = json.load(f)

    # what type of radar scan
    if "SUR" in radarfile:
        scan_type = "PPI"
    elif "RHI" in radarfile:
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
