import matplotlib
import polars as pl
import numpy as np
import metpy.calc as mpcalc
from metpy.plots import SkewT, Hodograph
from metpy.units import units
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

matplotlib.use("agg")


def plot_skewt(df, radiosonde_date, radiosonde_time, radiosonde_location, outdir):
    hgt = list(df["HeightMSL"]) * units.m
    p = list(df["P"]) * units.hPa
    T = list(df["Temp"]) * units.degC
    Td = list(df["Dewp"]) * units.degC
    wind_speed = list(df["Speed"] * 1.9438) * units.knots
    wind_dir = list(df["Dir"]) * units.degrees
    u, v = mpcalc.wind_components(wind_speed, wind_dir)

    fig = plt.figure(figsize=(9, 9))
    skew = SkewT(fig)
    skew.plot(p, T, "r")
    skew.plot(p, Td, "g")
    skew.plot_barbs(p[::80], u[::80], v[::80])

    skew.ax.set_xlabel(f"Temperature ({T.units:~P})")
    skew.ax.set_ylabel("Pressure (mb)")

    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    skew.plot(lcl_pressure, lcl_temperature, "ko", markerfacecolor="black")

    prof = mpcalc.parcel_profile(p, T[0], Td[0]).to("degC")
    skew.plot(p, prof, "k", linewidth=2)

    skew.shade_cin(p, T, prof, Td)
    skew.shade_cape(p, T, prof)

    skew.ax.set_xlim(-50, 30)
    skew.ax.set_ylim(max(p), 100)

    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    mixing_lines_pressure = units.Quantity(
        np.linspace(400, max(skew.ax.get_ylim())), "mbar"
    )
    skew.plot_mixing_lines(pressure=mixing_lines_pressure)

    ax_hod = inset_axes(skew.ax, "40%", "40%", loc=1)
    h = Hodograph(ax_hod, component_range=40.0)
    h.add_grid(increment=10)
    h.plot_colormapped(u, v, hgt)

    skew.ax.set_title(f"{radiosonde_location} - {radiosonde_date} {radiosonde_time}Z")

    if not os.path.exists(f"{outdir}/{radiosonde_date}"):
        os.mkdir(f"{outdir}/{radiosonde_date}")
    plt.savefig(
        f"{outdir}/{radiosonde_date}/skewt_{radiosonde_date}T{radiosonde_time}.png"
    )

    plt.close()


def main(infile, outdir):
    file_name = sys.argv[1]
    outdir = sys.argv[2]

    radiosonde_date = file_name.split("/")[-1].split("_")[1]
    radiosonde_time = file_name.split("/")[-1].split("_")[2]

    with open(file_name, encoding="charmap") as f:
        line_skip = 0
        units_line_next = False
        while True:
            data = f.readline()
            if units_line_next:
                # data_units = [i.strip() for i in data.split("\t")]
                break
            elif data.startswith("Elapsed time"):
                units_line_next = True
            elif "Station name" in data:
                radiosonde_location = "_".join(
                    data.split("Station name")[1].strip().split(" ")
                )
                line_skip += 1
            else:
                line_skip += 1

    df = pl.read_csv(
        file_name,
        encoding="charmap",
        skip_rows=line_skip,
        separator="\t",
        skip_rows_after_header=1,
        ignore_errors=True,
    )
    plot_skewt(df, radiosonde_date, radiosonde_time, radiosonde_location, outdir)


if __name__ == "__main__":
    import sys

    file_name = sys.argv[1]
    outdir = sys.argv[2]
    main(file_name, outdir)
