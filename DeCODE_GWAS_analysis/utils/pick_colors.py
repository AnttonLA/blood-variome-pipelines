import argparse
from matplotlib import cm

"""
This script generates a list of hexadecimal color values using the matplotlib colormap. It takes a single argument, the
number of colors to generate. It prints the list of hexadecimal color values as a comma-separated list, with no spaces.

Usage:
    python pick_colors.py 10
"""

# Create an argument parser object
parser = argparse.ArgumentParser()

# Add an argument for the number of colors to generate
parser.add_argument('num_colors', type=int, help='the number of colors to generate')

# Parse the command line arguments
args = parser.parse_args()

# Generate the color palette using the matplotlib colormap
colors = cm.get_cmap('twilight', args.num_colors)


def rgba_to_hex(rgba: tuple) -> str:
    # Convert the RGBA values to integers in the range 0-255
    r, g, b, a = [int(x * 255) for x in rgba]

    # Return the hexadecimal color value in the format #RRGGBB
    return f"#{r:02X}{g:02X}{b:02X}"


# Print the hexadecimal values as a comma-separated list, with no spaces
print(",".join([rgba_to_hex(colors(i)) for i in range(args.num_colors)]))
