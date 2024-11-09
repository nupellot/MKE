# Use x11 terminal for interactive display
set terminal wxt size 800,600 enhanced font 'Arial,10'

# Graph title
set title 'Comparison of Analytical Solution and FEM Solution'

# Axis labels
set xlabel 'x' # X-axis label
set ylabel 'u(x)' # Y-axis label

# Graph legend
set key top left

# Grid on the graph
set grid

# Style settings for lines
set style line 1 lt 1 lw 2 pt 7 ps 1.5 lc rgb 'blue'   # Style for analytical solution (blue dots)
set style line 2 lt 1 lw 2 pt 5 ps 1.5 lc rgb 'red'    # Style for FEM solution (red dots)

# Check for correct data reading
print "Reading data from 'Real.txt' and 'MKE.txt'"

# Plot the graph
plot 'Real.txt' using 1:2 with linespoints linestyle 1 title "Analitical Solution", \
     'MKE.txt' using 1:2 with linespoints linestyle 2 title 'FEM', \
     'AbsolutelyReal.txt' using 1:2 with linespoints title "AbsolutelyReal Solution"

# Wait for user input before closing the window
pause -1 "Press Enter to close the window"
