"""ani.py - A set of functions to animate particular results.

python src/visualisation/animate.py
"""
from typing import Callable
import numpy as np
import xarray as xr
from tqdm import tqdm
import matplotlib.pyplot as plt
import imageio
from src.plot_settings import ps_defaults  # label_subplots,
from src.utils import timeit
from src.read_in_data import ReadData


@timeit
def animate(
    video_path: str = "gifs/ani_model.mp4",
) -> None:
    """This function animates the inputs, labels, predictions.

    Based on code originally from Tom Anderson: tomand@bas.ac.uk.
    """
    ps_defaults(use_tex=False, dpi=200)

    read = ReadData()
    ds_uv, _ = read.get_time_datasets()

    def balance(tup):
        return min(tup[0], -tup[1]), max(-tup[0], tup[1])

    def mi(a):
        return np.percentile(a, 2)

    def ma(a):
        return np.percentile(a, 98)

    umin, umax = balance((mi(ds_uv.ucon.values), ma(ds_uv.vcon.values)))
    vmin, vmax = balance((mi(ds_uv.vcon.values), ma(ds_uv.vcon.values)))
    wmin, wmax = balance((mi(ds_uv.wcon.values), ma(ds_uv.wcon.values)))

    cmap = "seismic"

    def gen_frame_func() -> Callable:
        """Create imageio frame function for xarray.DataArray visualisation.

        Returns:
            make_frame (Callable): function to create each frame.

        """

        def make_frame(index: int) -> np.array:
            """Make an individual frame of the animation.

            Args:
                index (int): The time index.

            Returns:
                image (np.array): np.frombuffer output that can be fed into imageio
            """
            time = index
            fig, axs = plt.subplots(3, 1, sharex=True)
            ds_uv.ucon.isel(T=time, x=slice(0, 100)).plot(
                ax=axs[0],
                vmin=umin,
                vmax=umax,
                cmap=cmap,
                cbar_kwargs={"label": "$u$ [m s$^{-1}$]"},
            )
            axs[0].set_xlabel("")
            axs[0].set_title("Day " + "{:.0f}".format(ds_uv.coords["T"].values[time]))
            ds_uv.vcon.isel(T=time, x=slice(0, 100)).plot(
                ax=axs[1],
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
                cbar_kwargs={"label": "$v$ [m s$^{-1}$]"},
            )
            axs[1].set_xlabel("")
            axs[1].set_title("")
            ds_uv.wcon.isel(T=time, x=slice(0, 100)).plot(
                ax=axs[2],
                vmin=wmin,
                vmax=wmax,
                cmap=cmap,
                cbar_kwargs={"label": "$w$ [m s$^{-1}$]"},
            )
            axs[2].set_title("")
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype="uint8")
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            plt.close()
            return image

        return make_frame

    video_indices = list(range(9))
    make_frame = gen_frame_func()
    imageio.mimsave(
        video_path,
        [make_frame(index) for index in tqdm(video_indices, desc=video_path)],
        fps=2,
    )
    print("Video " + video_path + " made.")


if __name__ == "__main__":
    # python src/visualisation/ani.py
    animate()
