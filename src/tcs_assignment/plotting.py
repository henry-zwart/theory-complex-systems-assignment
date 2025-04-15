from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import font_manager
import seaborn as sns


def configure_mpl():
    """Configure Matplotlib style."""
    FONT_SIZE_SMALL = 8
    FONT_SIZE_DEFAULT = 10
    FONT_PATH = Path("fonts/LibertinusSerif-Regular.otf")
    font_manager.fontManager.addfont(FONT_PATH)

    plt.rc("font", family="Libertinus Serif")
    plt.rc("font", weight="normal")  # controls default font
    plt.rc("mathtext", fontset="stix")
    plt.rc("font", size=FONT_SIZE_DEFAULT)  # controls default text sizes
    plt.rc("axes", titlesize=FONT_SIZE_DEFAULT)  # fontsize of the axes title
    plt.rc("axes", labelsize=FONT_SIZE_DEFAULT)  # fontsize of the x and y labels
    plt.rc("figure", labelsize=FONT_SIZE_DEFAULT)
    plt.rc("figure", dpi=600)

    sns.set_context(
        "paper",
        rc={
            "axes.linewidth": 0.5,
            "axes.labelsize": FONT_SIZE_DEFAULT,
            "axes.titlesize": FONT_SIZE_DEFAULT,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "ytick.minor.width": 0.4,
            "xtick.labelsize": FONT_SIZE_SMALL,
            "ytick.labelsize": FONT_SIZE_SMALL,
        },
    )
