#!/bin/env python
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import sys

def plot_data_cube(cubedata):
    # Initialize the figure
    fig, ax = plt.subplots()
    fig.patch.set_facecolor('lightgray')  # Set the background color of the figure
    plt.subplots_adjust(bottom=0.25)  # Adjust the position of the slider
    ax_slider = plt.axes([0.1, 0.1, 0.8,.05])

    # Create a slider for selecting the frame
    frame_slider = Slider(ax_slider, 'Frame', 0, cubedata.shape[0] - 1, valinit=0, valstep=1)

    # Function to update the plot
    def update(val):
        frame = int(frame_slider.val)
        ax.clear()  # Clear the previous frame
        ax.imshow(cubedata[frame, :, :], cmap='viridis') # , aspect='auto'
        ax.set_title(f'Frame {frame}')
        plt.draw()

    frame_slider.on_changed(update)  # Connect the slider to the update function

    # Initial plot
    update(0)


    plt.show()

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot a 3D data cube")
    # pass path to .npy file "cubes_ugc89/tr_images/u89cl_054.npy"
    parser.add_argument("input_path", metavar="input_npy_path", type=str, help="Path to the input .npy file")
    
    args = parser.parse_args()

    bindata_path = args.input_path
    cubedata = np.load(bindata_path)
    cubedata = np.where(cubedata<0, np.zeros_like(cubedata), cubedata)

    plot_data_cube(cubedata)

if __name__ == "__main__":
    main()
