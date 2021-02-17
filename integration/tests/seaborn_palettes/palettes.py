"""Show seaborn palettes and simulate colorblind versions."""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from colorspacious import cspace_converter


def show_palettes():

    f = plt.figure(figsize=(6, 6))

    ax_locs = dict(
        deep=(.4, .4),
        bright=(.8, .8),
        muted=(.49, .71),
        dark=(.8, .2),
        pastel=(.2, .8),
        colorblind=(.71, .49),
    )

    s = .35

    for pal, (x, y) in ax_locs.items():
        ax = f.add_axes([x - s / 2, y - s / 2, s, s])
        ax.pie(np.ones(10),
               colors=sns.color_palette(pal, 10),
               counterclock=False, startangle=180,
               wedgeprops=dict(linewidth=1, edgecolor="w"))
        f.text(x, y, pal, ha="center", va="center", size=14,
               bbox=dict(facecolor="white", alpha=0.85,
                         boxstyle="round,pad=0.2"))

    f.text(.1, .05, "Saturation", size=18, ha="left", va="center",
           bbox=dict(facecolor="white", edgecolor="w"))
    f.text(.05, .1, "Luminance", size=18, ha="center", va="bottom",
           rotation=90, bbox=dict(facecolor="white", edgecolor="w"))

    ax = f.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.arrow(.15, .05, .4, 0, width=.002, head_width=.015, color="k")
    ax.arrow(.05, .15, 0, .4, width=.002, head_width=.015, color="k")

    return f


if __name__ == "__main__":

    f = show_palettes()
    f.suptitle("Original palettes", size=14)
    f.savefig("palettes.png", dpi=100)
    plt.close(f)

    converters = {}

    _deuter50_space = {"name": "sRGB1+CVD",
                       "cvd_type": "deuteranomaly",
                       "severity": 50}
    converters["deuter50"] = cspace_converter(_deuter50_space, "sRGB1")
    _deuter100_space = {"name": "sRGB1+CVD",
                        "cvd_type": "deuteranomaly",
                        "severity": 100}
    converters["deuter100"] = cspace_converter(_deuter100_space, "sRGB1")

    _prot50_space = {"name": "sRGB1+CVD",
                     "cvd_type": "protanomaly",
                     "severity": 50}
    converters["prot50"] = cspace_converter(_prot50_space, "sRGB1")

    _prot100_space = {"name": "sRGB1+CVD",
                      "cvd_type": "protanomaly",
                      "severity": 100}
    converters["prot100"] = cspace_converter(_prot100_space, "sRGB1")

    for cbt, converter in converters.items():
        print(cbt, converter)

        f = show_palettes()
        for ax in f.axes:
            for patch in ax.patches:
                fc = patch.get_facecolor()[:3]
                print(fc, converter(fc))
                patch.set_facecolor(np.clip(converter(fc), 0, 1))

        f.suptitle(cbt.capitalize(), size=14)
        f.savefig(f"palettes_{cbt}.png", dpi=100)
        plt.close(f)