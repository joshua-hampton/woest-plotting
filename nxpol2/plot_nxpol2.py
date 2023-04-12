import pyart
import cartopy.crs as ccrs

# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as lines

# from matplotlib.colors import ListedColormap
import json
import sys
from utils import wescon_kml_grid

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
    num_colours = var_scales[radar_name][var]["num_colours"]
    colourmap = var_scales[radar_name][var]["colourmap"]
    print(num_colours, colourmap)

    # radar display
    display = pyart.graph.RadarMapDisplay(radar)

    # make plots for each sweep
    for sweep in range(radar.nsweeps):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        display.plot_ppi_map(
            var,
            ax=ax,
            resolution="10m",
            vmin=vmin,
            vmax=vmax,
            sweep=sweep,
            # colorbar_orient="horizontal",
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
        sweep_elevation = radar.elevation["data"][
            sum(radar.rays_per_sweep["data"][0 : sweep + 1]) - 1
        ]
        plt.savefig(
            f"{outdir}/{radar_name}_{radar.metadata['start_datetime'].replace(('-','').replace(':',''))}_{var}_ele{sweep_elevation}.png"
        )
        plt.close()


def make_rhi_plots(radar, radar_name, outdir, var, var_scales):
    # get variable info
    vmin = var_scales[radar_name][var]["min"]
    vmax = var_scales[radar_name][var]["max"]
    num_colours = var_scales[radar_name][var]["num_colours"]
    colourmap = var_scales[radar_name][var]["colourmap"]
    return vmin, vmax, num_colours, colourmap


def main(
    radarfile, outdir, variables=variables, radar_name="ncas-mobile-x-band-radar-2"
):
    radar = pyart.io.read(radarfile)

    # read in stuff from utils
    outer_lines, h_lines, v_lines = wescon_kml_grid.read_kml("../WesConGrid.kml")
    with open("utils/var_scales.json") as f:
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
