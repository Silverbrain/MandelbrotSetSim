# Set the output file type and name
set terminal pngcairo enhanced size 12000, 8000
set output 'mandelbrot.png'

# Remove axes, labels, and title
unset border
unset xtics
unset ytics
unset xlabel
unset ylabel
unset title

# Configure the plot appearance
set size ratio -1       # Ensures correct aspect ratio
set palette defined (0 "black", 1 "blue", 2 "red", 3 "yellow", 4 "white")  # Custom palette
unset colorbox          # Remove color legend (optional)

# Plot the Mandelbrot data as a heatmap
plot 'mandelbrot_sim.dat' using 1:2:3 with image